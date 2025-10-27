
struct TPT1msaFunctional{P} <: AbstractFunctional
    weight_fft_K :: Array{Float64, P}
    model        :: TPT1MSAFreeEnergy
end

"""
    construct_functional(::Type{TPT1msaFunctional}, molsys::MolecularSystem, bulk::BulkState, geometry::CartesianCoord) -> TPT1msaFunctional

Constructs a `TPT1msaFunctional` object for electrostatic correlation using TPT1 + MSA model.

# Arguments
- `molsys`: Molecular system containing monomer diameters.
- `bulk`: Bulk state with free energy model instances.
- `geometry`: Defines spatial discretization (`NP`, `bin_width`).

# Returns
- A fully initialized `TPT1msaFunctional` object.
"""
function construct_functional(::Type{TPT1msaFunctional},
                              molsys::MolecularSystem,
                              bulk::BulkState,
                              geometry::CartesianCoord)
    @unpack NP, bin_width, mirrored, offset = geometry
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
        @inbounds for K_star in Rsys_star
            weight_fft_K[K_star, j] = weight_fft[K_star] / (4π * bval^2)
        end
    end

    model = only(filter(e -> isa(e, TPT1MSAFreeEnergy), bulk.fe_model))

    return TPT1msaFunctional{ndims(weight_fft_K)}(weight_fft_K, model)
end



function compute_bead_weighted_densities(bulk_system::IsingLST,
                                         fields::SpatialFields,
                                         fe_term::TPT1msaFunctional)
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


function compute_bond_weighted_densities(bulk_system::IsingLST,
                                              fields::SpatialFields,
                                              fe_term::TPT1msaFunctional)
    @unpack rho_bonds_hat, f_hat, plan_backward = fields.fourier
    @unpack weight_fft_K = fe_term
    Rsys_star = CartesianIndices(f_hat)

    wt_bonds_K = similar(rho_bonds_hat, Float64)

    # Bond weighted densities
    for (u, bond) in enumerate(bulk_system.molsys.properties.bond_types)
        for (idx, i) in enumerate(bond)
            for K_star in Rsys_star
                f_hat[K_star] = rho_bonds_hat[K_star, idx, u] * weight_fft_K[K_star, i]
            end

            f_hat = plan_backward * f_hat

            for K_star in Rsys_star
                wt_bonds_K[K_star, idx, u] = abs(real(f_hat[K_star]))
            end
        end
    end

    return wt_bonds_K
end


function compute_pointwise_free_energy(fe_term::TPT1msaFunctional,
                                       wt_rho_K::AbstractArray,
                                       wt_bonds_K::AbstractArray,
                                       bulk_system, geometry)
    @unpack NP, mirrored, offset = geometry

    @unpack wt_rho, wt_bonds = fe_term.model
    @unpack bond_types = bulk_system.molsys.properties
    @unpack diameters, valences = bulk_system.molsys.properties.monomers
    @unpack bjerrum_length = bulk_system.molsys.properties.system

    fe = 0.0
    for K in CartesianIndices(NP)
        K_star = to_star_index(K, offset, mirrored)

        @views @. wt_rho   = wt_rho_K[K_star, :]
        @views @. wt_bonds = wt_bonds_K[K_star, :, :]

        fe += free_energy_TPT1msa(wt_bonds, bond_types, wt_rho,
                                  diameters, valences, bjerrum_length)
    end

    return fe
end


function compute_pointwise_derivative(fe_term::TPT1msaFunctional,
                                      wt_rho_K::AbstractArray,
                                      wt_bonds_K::AbstractArray,
                                      bulk_system)

    dPhi_K = similar(wt_rho_K)
    lng_K  = similar(wt_bonds_K)

    @unpack wt_rho, wt_bonds, dPhi, lng = fe_term.model
    @unpack bond_types = bulk_system.molsys.properties
    @unpack diameters, valences = bulk_system.molsys.properties.monomers
    @unpack bjerrum_length = bulk_system.molsys.properties.system

    for K_star in CartesianIndices(axes(wt_rho_K)[1:end-1])
        @views @. wt_rho   = wt_rho_K[K_star, :]
        @views @. wt_bonds = wt_bonds_K[K_star, :, :]

        compute_TPT1msa_deriv!(dPhi, lng, wt_bonds, bond_types,
                               wt_rho, diameters, valences, bjerrum_length)

        @views @. dPhi_K[K_star, :] = dPhi
        for u in eachindex(lng)
            @views @. lng_K[K_star, :, u] = lng[u]
        end
    end

    return dPhi_K, lng_K
end


function accumulate_mu_ex!(bulk_system::IsingLST,
                           fields::SpatialFields,
                           fe_term::TPT1msaFunctional,
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


function accumulate_lng!(bulk_system::IsingLST,
                         fields::SpatialFields,
                         fe_term::TPT1msaFunctional,
                         lng_K::AbstractArray)

    @unpack f_hat, plan_forward, lng_hat = fields.fourier
    @unpack weight_fft_K = fe_term
    @unpack bond_types = bulk_system.molsys.properties

    Rsys_star = CartesianIndices(f_hat)

    for (u, bond) in enumerate(bond_types)
        for (idx, i) in enumerate(bond)
            for K_star in Rsys_star
                f_hat[K_star] = lng_K[K_star, idx, u]
            end

            f_hat = plan_forward * f_hat

            for K_star in Rsys_star
                lng_hat[K_star, idx, u] += f_hat[K_star] * weight_fft_K[K_star, i]
            end
        end
    end
end


function eval_fe_functional(bulk_system::IsingLST,
                            geometry::CoordSystem,
                            fields::SpatialFields,
                            fe_term::TPT1msaFunctional)
    wt_rho_K = compute_bead_weighted_densities(bulk_system, fields, fe_term)
    wt_bonds_K = compute_bond_weighted_densities(bulk_system, fields, fe_term)

    return compute_pointwise_free_energy(fe_term, wt_rho_K, wt_bonds_K, bulk_system, geometry)
end


function eval_mu_functional!(bulk_system::IsingLST,
                             fields::SpatialFields,
                             fe_term::TPT1msaFunctional)
    wt_rho_K = compute_bead_weighted_densities(bulk_system, fields, fe_term)
    wt_bonds_K = compute_bond_weighted_densities(bulk_system, fields, fe_term)

    dPhi_K, lng_K = compute_pointwise_derivative(fe_term, wt_rho_K, wt_bonds_K, bulk_system)

    accumulate_mu_ex!(bulk_system, fields, fe_term, dPhi_K)

    accumulate_lng!(bulk_system, fields, fe_term, lng_K)
end
