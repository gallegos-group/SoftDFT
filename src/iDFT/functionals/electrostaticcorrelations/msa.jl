struct msaFunctional{O} <: AbstractFunctional
    weight_fft_K :: Array{Float64, O}
    model   :: MSAFreeEnergy
end


"""
    construct_functional(::Type{msaFunctional}, molsys::MolecularSystem, bulk::BulkState, geometry::CartesianCoord) -> msaFunctional

Constructs an `msaFunctional` object using the system’s geometry and model parameters.

# Arguments
- `molsys`: Molecular system containing monomer properties (e.g., diameters).
- `bulk`: Bulk state that holds initialized free energy models.
- `geometry`: Defines the spatial discretization of the system (`NP`, `bin_width`).

# Returns
- A fully initialized `msaFunctional` object.
"""
function construct_functional(::Type{msaFunctional}, molsys::MolecularSystem, bulk::BulkState, geometry::CartesianCoord)
    @unpack NP, bin_width, features = geometry
    @unpack mirrored, offset = features
    @unpack diameters = molsys.properties.monomers

    # work on the calculation domain
    NP_star   = compute_full_domain(NP, mirrored, offset)
    dims_star = (NP_star..., length(diameters))

    weight_fft_K = zeros(Float64, dims_star)
    weight_fft   = zeros(Float64, NP_star)
    Rsys_star    = CartesianIndices(NP_star)

    for (j, σ) in enumerate(diameters)
        bval = σ / 2.0
        delta_func_fft!(weight_fft, bin_width, bval)
        for K_star in Rsys_star
            weight_fft_K[K_star, j] = weight_fft[K_star] / (4π * bval^2)
        end
    end

    model = only(filter(e -> isa(e, MSAFreeEnergy), bulk.fe_model))

    return msaFunctional{ndims(weight_fft_K)}(weight_fft_K, model)
end


function compute_bead_weighted_densities(bulk_system::IsingLST,
                                         fields::SpatialFields,
                                         fe_term::msaFunctional)
    @unpack rho_beads_hat, f_hat, plan_backward = fields.fourier
    @unpack weight_fft_K = fe_term
    Rsys_star = CartesianIndices(f_hat)

    wt_rho_K = similar(rho_beads_hat, Float64)

    for j in axes(rho_beads_hat, ndims(rho_beads_hat))
        for K_star in Rsys_star
            f_hat[K_star] = rho_beads_hat[K_star, j] * weight_fft_K[K_star, j]
        end

        f_hat = plan_backward * f_hat

        for K_star in Rsys_star
            wt_rho_K[K_star, j] = abs(real(f_hat[K_star]))
        end
    end

    return wt_rho_K
end


function compute_pointwise_free_energy(fe_term::msaFunctional,
                                       wt_rho_K::AbstractArray,
                                       bulk_system, geometry)
    @unpack NP, features = geometry
    @unpack mirrored, offset = features

    @unpack wt_rho = fe_term.model
    @unpack diameters, valences = bulk_system.molsys.properties.monomers
    @unpack bjerrum_length = bulk_system.molsys.properties.system

    fe = 0.0
    for K in CartesianIndices(NP)
        K_star = to_star_index(K, offset, mirrored)
        @views @. wt_rho = wt_rho_K[K_star, :]

        fe += free_energy_msa(wt_rho, diameters, valences, bjerrum_length)
    end

    return fe
end


function compute_pointwise_derivative(fe_term::msaFunctional,
                                      wt_rho_K::AbstractArray,
                                      bulk_system)

    dPhi_K = similar(wt_rho_K)

    @unpack dPhi, wt_rho = fe_term.model
    @unpack diameters, valences = bulk_system.molsys.properties.monomers
    @unpack bjerrum_length = bulk_system.molsys.properties.system

    for K_star in CartesianIndices(axes(wt_rho_K)[1:end-1])
        @views @. wt_rho = wt_rho_K[K_star, :]

        compute_msa_deriv!(dPhi, wt_rho, diameters, valences, bjerrum_length)

        @views @. dPhi_K[K_star, :] = dPhi
    end

    return dPhi_K
end


function accumulate_mu_ex!(bulk_system::IsingLST,
                           fields::SpatialFields,
                           fe_term::msaFunctional,
                           dPhi_K::AbstractArray)

    @unpack f_hat, plan_forward, mu_ex_hat = fields.fourier
    @unpack weight_fft_K = fe_term

    Rsys_star = CartesianIndices(f_hat)

    for j in axes(dPhi_K, ndims(dPhi_K))
        for K_star in Rsys_star
            f_hat[K_star] = dPhi_K[K_star, j]
        end

        f_hat = plan_forward * f_hat

        for K_star in Rsys_star
            mu_ex_hat[K_star, j] += f_hat[K_star] * weight_fft_K[K_star, j]
        end
    end
end


function eval_fe_functional(bulk_system::IsingLST,
                            geometry::CoordSystem,
                            fields::SpatialFields,
                            fe_term::msaFunctional)

    wt_rho_K = compute_bead_weighted_densities(bulk_system, fields, fe_term)

    return compute_pointwise_free_energy(fe_term, wt_rho_K, bulk_system, geometry)
end

function eval_mu_functional!(bulk_system::IsingLST,
                             fields::SpatialFields,
                             fe_term::msaFunctional)

    wt_rho_K = compute_bead_weighted_densities(bulk_system, fields, fe_term)

    dPhi_K   = compute_pointwise_derivative(fe_term, wt_rho_K, bulk_system)

    accumulate_mu_ex!(bulk_system, fields, fe_term, dPhi_K)
end
