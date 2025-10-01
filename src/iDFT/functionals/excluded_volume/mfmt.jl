struct mfmtFunctional{P} <: AbstractFunctional
    weight_fft_K :: Array{Float64, P}
    model        :: MFMTFreeEnergy
end

"""
    construct_functional(::Type{mfmtFunctional}, molsys::MolecularSystem, bulk::BulkState, geometry::CartesianCoord) -> mfmtFunctional

Constructs an `mfmtFunctional` representing the modified fundamental measure theory (MFMT) contribution.

# Arguments
- `molsys`: Molecular system with monomer properties (e.g., diameters).
- `bulk`: Bulk state containing the initialized MFMT model.
- `geometry`: Spatial discretization (`NP`, `bin_width`).

# Returns
- A fully initialized `mfmtFunctional` object.
"""
function construct_functional(::Type{mfmtFunctional}, molsys::MolecularSystem, bulk::BulkState, geometry::CartesianCoord)
    @unpack NP, bin_width, features = geometry
    @unpack mirrored, offset = features
    @unpack diameters = molsys.properties.monomers

    # work on the calculation domain
    NP_star   = compute_full_domain(NP, mirrored, offset)
    dims_star = (NP_star..., 5, length(diameters))

    weight_fft_K = zeros(Float64, dims_star)
    weight_fft   = zeros(Float64, NP_star)
    Rsys_star    = CartesianIndices(NP_star)

    # Scalar weight: n2 (delta)
    for (j, σ) in enumerate(diameters)
        bval = σ / 2.0
        delta_func_fft!(weight_fft, bin_width, bval)
        for K_star in Rsys_star
            weight_fft_K[K_star, 1, j] = weight_fft[K_star]
        end
    end

    # Scalar weight: n3 (step)
    for (j, σ) in enumerate(diameters)
        bval = σ / 2.0
        step_func_fft!(weight_fft, bin_width, bval)
        for K_star in Rsys_star
            weight_fft_K[K_star, 2, j] = weight_fft[K_star]
        end
    end

    # Vector weights: nV2x, nV2y, nV2z (directional delta)
    for (j, σ) in enumerate(diameters)
        bval = σ / 2.0
        for v in eachindex(NP)
            vectored_delta_func_fft!(weight_fft, v, bin_width, bval)
            for K_star in Rsys_star
                weight_fft_K[K_star, v + 2, j] = weight_fft[K_star]
            end
        end
    end

    model = only(filter(e -> isa(e, MFMTFreeEnergy), bulk.fe_model))
    return mfmtFunctional{ndims(weight_fft_K)}(weight_fft_K, model)
end


function compute_bead_weighted_densities(bulk_system::IsingLST,
                                         fields::SpatialFields,
                                         fe_term::mfmtFunctional)

    @unpack rho_beads_hat, plan_backward, f_hat = fields.fourier
    @unpack weight_fft_K = fe_term

    NB = size(rho_beads_hat, ndims(rho_beads_hat))
    nai_K = similar(f_hat, Float64, size(f_hat)..., 5, NB)
    Rsys_star = CartesianIndices(f_hat)

    for j in axes(rho_beads_hat, ndims(rho_beads_hat))
        # scalar (m=1,2)
        for m in 1:2
            @. f_hat = 0.0
            @inbounds for K_star in Rsys_star
                f_hat[K_star] = rho_beads_hat[K_star, j] * weight_fft_K[K_star, m, j]
            end
            f_hat = plan_backward * f_hat
            @inbounds for K_star in Rsys_star
                nai_K[K_star, m, j] = abs(real(f_hat[K_star]))
            end
        end

        # vector (m=3:5)
        for m in 3:5
            @. f_hat = 0.0
            @inbounds for K_star in Rsys_star
                f_hat[K_star] = rho_beads_hat[K_star, j] * (im * weight_fft_K[K_star, m, j])
            end
            f_hat = plan_backward * f_hat
            @inbounds for K_star in Rsys_star
                nai_K[K_star, m, j] = real(f_hat[K_star])
            end
        end
    end

    return nai_K
end


function pointwise_weighted_densities!(nai_K, K_star, diameters, nV1, nV2)
    n0 = 0.0
    n1 = 0.0
    n2 = 0.0
    n3 = 0.0
    @. nV1 = 0.0
    @. nV2 = 0.0

    for (j, σ) in enumerate(diameters)
        n0 += nai_K[K_star, 1, j] / (π*σ^2)
        n1 += nai_K[K_star, 1, j] / (2π*σ)
        n2 += nai_K[K_star, 1, j]
        n3 += nai_K[K_star, 2, j]
        for v in eachindex(nV1)
            nV1[v] += nai_K[K_star, v+2, j] / (2π*σ)
            nV2[v] += nai_K[K_star, v+2, j]
        end
    end

    # cap at 0.95
    if n3 > 0.95
        fac = 0.95 / n3
        n0 *= fac; n1 *= fac; n2 *= fac; n3 *= fac
        @. nV1 *= fac
        @. nV2 *= fac
    end

    return n0, n1, n2, n3
end


function compute_pointwise_free_energy(fe_term::mfmtFunctional,
                                       nai_K::AbstractArray,
                                       bulk_system, geometry)
    @unpack NP, features = geometry
    @unpack mirrored, offset = features

    @unpack nV1, nV2 = fe_term.model
    @unpack diameters = bulk_system.molsys.properties.monomers

    fe = 0.0
    for K in CartesianIndices(NP)
        K_star = to_star_index(K, offset, mirrored)
        n0, n1, n2, n3 = pointwise_weighted_densities!(nai_K, K_star, diameters, nV1, nV2)

        fe += free_energy_mfmt(n0, n1, n2, n3, nV1, nV2)
    end

    return fe
end


function eval_fe_functional(bulk_system :: IsingLST, geometry :: CoordSystem, fields :: SpatialFields, fe_term :: mfmtFunctional)

    @unpack rho_beads_hat, plan_backward, f_hat = fields.fourier
    @unpack weight_fft_K = fe_term

    # Calculate weighted densities
    nai_K = compute_bead_weighted_densities(bulk_system, fields, fe_term)

    return compute_pointwise_free_energy(fe_term, nai_K, bulk_system, geometry)
end


function compute_pointwise_derivative(fe_term::mfmtFunctional,
                                      nai_K::AbstractArray,
                                      bulk_system)

    dPhi_K = similar(nai_K)  # allocate output

    @unpack dPhi, nV1, nV2 = fe_term.model
    @unpack diameters = bulk_system.molsys.properties.monomers

    for K_star in CartesianIndices(axes(dPhi_K)[1:end-2])
        n0, n1, n2, n3 = pointwise_weighted_densities!(nai_K, K_star, diameters, nV1, nV2)

        # derivative evaluation
        if n3 > 1e-10
            compute_mfmt_deriv!(dPhi, n0, n1, n2, n3, nV1, nV2)
        else
            @. dPhi = 0.0
        end

        # write back into dPhi_K
        @views @. dPhi_K[K_star, 1, :] = dPhi[1]/(π*diameters^2) +
                                         dPhi[2]/(2π*diameters) +
                                         dPhi[3]
        @views @. dPhi_K[K_star, 2, :] = dPhi[4]
        for v in eachindex(nV1)
            @views @. dPhi_K[K_star, v+2, :] = dPhi[v+4]/(2π*diameters) + dPhi[v+7]
        end
    end

    return dPhi_K
end


function accumulate_mu_ex!(bulk_system::IsingLST,
                           fields::SpatialFields,
                           fe_term::mfmtFunctional,
                           dPhi_K)

    @unpack f_hat, plan_forward, mu_ex_hat = fields.fourier
    @unpack weight_fft_K = fe_term

    Rsys_star = CartesianIndices(f_hat)

    for j in axes(dPhi_K, ndims(dPhi_K))   # last dim = bead index
        # scalars (m = 1,2)
        for m = 1:2
            @inbounds for K_star in Rsys_star
                f_hat[K_star] = dPhi_K[K_star, m, j]
            end
            f_hat = plan_forward * f_hat
            @inbounds for K_star in Rsys_star
                mu_ex_hat[K_star, j] += f_hat[K_star] * weight_fft_K[K_star, m, j]
            end
        end

        # vectors (m = 3:5)
        for m = 3:5
            @inbounds for K_star in Rsys_star
                f_hat[K_star] = dPhi_K[K_star, m, j]
            end
            f_hat = plan_forward * f_hat
            @inbounds for K_star in Rsys_star
                mu_ex_hat[K_star, j] += f_hat[K_star] * (-im * weight_fft_K[K_star, m, j])
            end
        end
    end

    return mu_ex_hat
end


function eval_mu_functional!(bulk_system :: IsingLST, fields :: SpatialFields, fe_term :: mfmtFunctional)

    # Calculate weighted densities
    nai_K = compute_bead_weighted_densities(bulk_system, fields, fe_term)

    # Calculate point-wise derivative
    dPhi_K = compute_pointwise_derivative(fe_term, nai_K, bulk_system)

    # Calculate excess chemical contribution by weighting point-wise derivative    
    accumulate_mu_ex!(bulk_system, fields, fe_term, dPhi_K)
end