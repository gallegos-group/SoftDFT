
struct TPT1mfmtFunctional{P} <: AbstractFunctional
    weight_fft_K :: Array{Float64, P}
    model        :: TPT1MFMTFreeEnergy
end

"""
    construct_functional(::Type{TPT1mfmtFunctional}, molsys::MolecularSystem, bulk::BulkState, geometry::CartesianCoord) -> TPT1mfmtFunctional

Constructs a `TPT1mfmtFunctional` object for the TPT1-MFMT free energy model.

# Arguments
- `molsys`: Molecular system including segment diameters.
- `bulk`: Bulk state containing the TPT1-MFMT model.
- `geometry`: Simulation grid geometry (`NP`, `bin_width`).

# Returns
- A fully initialized `TPT1mfmtFunctional` object.
"""
function construct_functional(::Type{TPT1mfmtFunctional}, molsys::MolecularSystem, bulk::BulkState, geometry::CartesianCoord)
    @unpack NP, bin_width, features = geometry
    @unpack mirrored, offset = features
    @unpack diameters = molsys.properties.monomers

    # work on the calculation domain
    NP_star   = compute_full_domain(NP, mirrored, offset)
    dims_star = (NP_star..., 5, length(diameters))

    weight_fft_K = zeros(Float64, dims_star)
    weight_fft   = zeros(Float64, NP_star)
    Rsys_star    = CartesianIndices(NP_star)

    # n2 term: scalar surface area
    for (j, σ) in enumerate(diameters)
        bval = σ / 2.0
        delta_func_fft!(weight_fft, bin_width, bval)
        for K_star in Rsys_star
            weight_fft_K[K_star, 1, j] = weight_fft[K_star]
        end
    end

    # n3 term: scalar volume
    for (j, σ) in enumerate(diameters)
        bval = σ / 2.0
        step_func_fft!(weight_fft, bin_width, bval)
        for K_star in Rsys_star
            weight_fft_K[K_star, 2, j] = weight_fft[K_star]
        end
    end

    # nV2 vector terms (x, y, z)
    for (j, σ) in enumerate(diameters)
        bval = σ / 2.0
        for v in eachindex(NP)
            bval = diameters[j] / 2.0
            vectored_delta_func_fft!(weight_fft, v, bin_width, bval)
            for K_star in Rsys_star
                weight_fft_K[K_star, v + 2, j] = weight_fft[K_star]
            end
        end
    end

    model = only(filter(e -> isa(e, TPT1MFMTFreeEnergy), bulk.fe_model))
    return TPT1mfmtFunctional{ndims(weight_fft_K)}(weight_fft_K, model)
end


function compute_bead_weighted_densities(bulk_system::IsingLST,
                                         fields::SpatialFields,
                                         fe_term::TPT1mfmtFunctional)

    @unpack rho_beads_hat, plan_backward, f_hat = fields.fourier
    @unpack weight_fft_K = fe_term

    NB = size(rho_beads_hat, ndims(rho_beads_hat))
    nai_K = similar(f_hat, Float64, size(f_hat)..., 5, NB)
    Rsys_star = CartesianIndices(f_hat)

    for j in axes(rho_beads_hat, ndims(rho_beads_hat))
        # scalar (m=1,2)
        for m in 1:2
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


function compute_bond_weighted_densities(bulk_system::IsingLST,
                                         fields::SpatialFields,
                                         fe_term::TPT1mfmtFunctional)

    @unpack rho_bonds_hat, plan_backward, f_hat = fields.fourier
    @unpack weight_fft_K = fe_term
    @unpack diameters = bulk_system.molsys.properties.monomers
    bond_types = bulk_system.molsys.properties.bond_types

    wt_bonds_K = similar(rho_bonds_hat, Float64)
    Rsys_star = CartesianIndices(f_hat)

    for (u, bond) in enumerate(bond_types)
        for (idx, i) in enumerate(bond)
            @inbounds for K_star in Rsys_star
                f_hat[K_star] = rho_bonds_hat[K_star, idx, u] * weight_fft_K[K_star, 1, i] / diameters[i]^2 / π
            end

            f_hat = plan_backward * f_hat

            @inbounds for K_star in Rsys_star
                wt_bonds_K[K_star, idx, u] = abs(real(f_hat[K_star]))
            end
        end
    end

    return wt_bonds_K
end

function pointwise_weighted_densities_i!(nai_K, K_star, diameters,
                                            n0i, n1i, n2i, n3i, nV1, nV2)
    # reset accumulators
    @. nV1 = 0.0
    @. nV2 = 0.0

    for (j, σ) in enumerate(diameters)
        n0i[j] = nai_K[K_star, 1, j] / (π*σ^2)
        n1i[j] = nai_K[K_star, 1, j] / (2π*σ)
        n2i[j] = nai_K[K_star, 1, j]
        n3i[j] = nai_K[K_star, 2, j]

        for v in eachindex(nV1)
            nV1[v] += nai_K[K_star, v+2, j] / (2π*σ)
            nV2[v] += nai_K[K_star, v+2, j]
        end
    end

    # cap at 0.95
    if sum(n3i) > 0.95
        fac = 0.95 / sum(n3i)
        @. n0i *= fac; @. n1i *= fac; @. n2i *= fac; @. n3i *= fac
        @. nV1 *= fac
        @. nV2 *= fac
    end
end


function compute_pointwise_free_energy(fe_term::TPT1mfmtFunctional,
                                       nai_K::AbstractArray,
                                       wt_bonds_K::AbstractArray,
                                       bulk_system, geometry)
    @unpack NP, features = geometry
    @unpack mirrored, offset = features

    @unpack wt_bonds = fe_term.model
    @unpack n0i, n1i, n2i, n3i, nV1, nV2 = fe_term.model
    @unpack diameters = bulk_system.molsys.properties.monomers
    @unpack bond_types = bulk_system.molsys.properties

    fe = 0.0
    for K in CartesianIndices(NP)
        K_star = to_star_index(K, offset, mirrored)
        @views @. wt_bonds = wt_bonds_K[K_star, :, :]

        pointwise_weighted_densities_i!(nai_K, K_star, diameters,
                                   n0i, n1i, n2i, n3i, nV1, nV2)

        fe += free_energy_TPT1mfmt(wt_bonds, bond_types, n0i, n1i, n2i, n3i, nV1, nV2, diameters)
    end

    return fe
end

function eval_fe_functional(bulk_system :: IsingLST, geometry :: CoordSystem, fields :: SpatialFields, fe_term :: TPT1mfmtFunctional)

    nai_K = compute_bead_weighted_densities(bulk_system, fields, fe_term)
    wt_bonds_K = compute_bond_weighted_densities(bulk_system, fields, fe_term)

    return compute_pointwise_free_energy(fe_term, nai_K, wt_bonds_K, bulk_system, geometry)
end


function compute_pointwise_derivative(fe_term::TPT1mfmtFunctional,
                                      nai_K::AbstractArray,
                                      wt_bonds_K::AbstractArray,
                                      bulk_system)

    dPhi_K = similar(nai_K)       # allocate output
    lng_K  = similar(wt_bonds_K)  # allocate output

    @unpack wt_bonds, dPhi, lng, n0i, n1i, n2i, n3i, nV1, nV2 = fe_term.model
    @unpack diameters = bulk_system.molsys.properties.monomers
    @unpack bond_types = bulk_system.molsys.properties

   for K_star in CartesianIndices(axes(dPhi_K)[1:end-2])
        @views @. wt_bonds = wt_bonds_K[K_star, :, :]

        pointwise_weighted_densities_i!(nai_K, K_star, diameters,
                                   n0i, n1i, n2i, n3i, nV1, nV2)

        # derivative evaluation
        if sum(n3i) > 1e-10
            compute_TPT1mfmt_deriv!(dPhi, lng, wt_bonds, bond_types,
                                    n0i, n1i, n2i, n3i, nV1, nV2, diameters)
        else
            @. dPhi = 0.0
            @. lng  = 0.0
        end

        # write back into dPhi_K
        @views @. dPhi_K[K_star, 1, :] = dPhi[1, :] / (π*diameters^2) +
                                          dPhi[2, :] / (2π*diameters) +
                                          dPhi[3, :]
        @views @. dPhi_K[K_star, 2, :] = dPhi[4, :]
        for v in 1:3
            @views @. dPhi_K[K_star, v+2, :] = dPhi[v+4, :] / (2π*diameters) +
                                                dPhi[v+7, :]
        end

        # write bond free energy logs
        for u in eachindex(lng)
            @views @. lng_K[K_star, :, u] = lng[u]
        end
    end

    return dPhi_K, lng_K
end


function accumulate_mu_ex!(bulk_system::IsingLST,
                           fields::SpatialFields,
                           fe_term::TPT1mfmtFunctional,
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


function accumulate_lng!(bulk_system::IsingLST,
                         fields::SpatialFields,
                         fe_term::TPT1mfmtFunctional,
                         lng_K::AbstractArray)

    @unpack f_hat, plan_forward, lng_hat = fields.fourier
    @unpack weight_fft_K = fe_term
    @unpack diameters = bulk_system.molsys.properties.monomers
    @unpack bond_types = bulk_system.molsys.properties

    Rsys_star = CartesianIndices(f_hat)

    for (u, bond) in enumerate(bond_types)
        for (idx, i) in enumerate(bond) 
            for K_star in Rsys_star
                f_hat[K_star] = lng_K[K_star, idx, u]
            end

            f_hat = plan_forward * f_hat

            for K_star in Rsys_star
                lng_hat[K_star, idx, u] += f_hat[K_star] * weight_fft_K[K_star, 1, i] / (π*diameters[i]^2)
            end
        end
    end

    return lng_hat
end


function eval_mu_functional!(bulk_system :: IsingLST, fields :: SpatialFields, fe_term :: TPT1mfmtFunctional)

    nai_K = compute_bead_weighted_densities(bulk_system, fields, fe_term)
    wt_bonds_K = compute_bond_weighted_densities(bulk_system, fields, fe_term)

    # Calculate point-wise derivative
    dPhi_K, lng_K = compute_pointwise_derivative(fe_term, nai_K, wt_bonds_K, bulk_system)
    
    # Calculate excess chemical contribution by weighting point-wise derivative    
    accumulate_mu_ex!(bulk_system, fields, fe_term, dPhi_K)

    # Bond contribution from ln g terms
    accumulate_lng!(bulk_system, fields, fe_term, lng_K)
end