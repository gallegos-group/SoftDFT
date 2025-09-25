"""
    struct coulFunctional{N} <: AbstractFunctional

Inhomogeneous DFT functional for the Coulomb contribution.

# Type Parameters
- `N`: Dimensionality of the simulation grid.

# Fields
- `weight_fft_K::Array{Float64,N}`: Fourier-space Coulomb kernel.
- `model::CoulombFreeEnergy`: Precomputed bulk free energy model.
"""
struct coulFunctional{N} <: AbstractFunctional
    weight_fft_K :: Array{Float64, N}
    model        :: CoulombFreeEnergy
end

"""
    construct_functional(::Type{coulFunctional}, molsys::MolecularSystem, bulk::BulkState, geometry::CartesianCoord) -> coulFunctional

Constructs a `coulFunctional` object using the provided molecular system, bulk state, and spatial geometry.

# Arguments
- `molsys`: Molecular system containing species, charges, and system parameters.
- `bulk`: Bulk state containing initialized free energy models and densities.
- `geometry`: Cartesian coordinate system specifying grid size and spacing.

# Returns
- A fully initialized `coulFunctional` object.
"""
function construct_functional(::Type{coulFunctional},
                              molsys::MolecularSystem,
                              bulk::BulkState,
                              geometry::CartesianCoord)
    @unpack NP, bin_width, features = geometry
    @unpack mirrored, offset = features
    @unpack bjerrum_length = molsys.properties.system

    # work on the calculation domain
    NP_star = compute_full_domain(NP, mirrored, offset)

    weight_fft_K = zeros(Float64, NP_star)
    coul_func_fft!(weight_fft_K, bin_width, bjerrum_length)

    model = only(filter(e -> isa(e, CoulombFreeEnergy), bulk.fe_model))

    return coulFunctional{ndims(weight_fft_K)}(weight_fft_K, model)
end




function eval_fe_functional(bulk_system :: IsingLST, geometry :: CoordSystem, fields :: SpatialFields, fe_term :: coulFunctional)
    @unpack features = geometry
    @unpack mirrored, offset = features

    @unpack weight_fft_K = fe_term
    @unpack surf_hat, rho_beads_hat, f_hat, mu_ex_hat, plan_backward = fields.fourier
    @unpack valences = bulk_system.molsys.properties.monomers
    PsiC = fields.excess.PsiC[1]
    Psi = fields.excess.Psi

    Rsys_star = CartesianIndices(f_hat)

    # Calculate excess chemical contribution by weighting point-wise derivative
    @. f_hat = surf_hat
    for j in axes(rho_beads_hat, ndims(rho_beads_hat))
        for K_star in Rsys_star
            f_hat[K_star] += rho_beads_hat[K_star, j] * valences[j]
        end

        # f_hat[1] += PsiC * valences[j]
    end

    @. f_hat *= weight_fft_K
    
    f_hat = plan_backward * f_hat

    for K in CartesianIndices(Psi)
        K_star = to_star_index(K, offset, mirrored)
        Psi[K] = real(f_hat[K_star])
    end
    
    @. f_hat = surf_hat
    for j in axes(rho_beads_hat, ndims(rho_beads_hat))
        for K_star in Rsys_star
            f_hat[K_star] += rho_beads_hat[K_star, j] * valences[j]
        end
    end
    
    f_hat = plan_backward * f_hat

    @. Psi += PsiC

    fe = 0.0
    for K in CartesianIndices(Psi)
        K_star = to_star_index(K, offset, mirrored)
        fe += real(f_hat[K_star]) * Psi[K]/ 2.0
    end
    # for K in CartesianIndices(Psi)
    #     K_star = to_star_index(K, offset, mirrored)
    #     fe += real(f_hat[K_star]) * Psi[K] / 2.0
    # end

    return fe
end


function eval_mu_functional!(bulk_system :: IsingLST, fields :: SpatialFields, fe_term :: coulFunctional)

    @unpack weight_fft_K = fe_term
    @unpack surf_hat, rho_beads_hat, f_hat, mu_ex_hat, plan_backward = fields.fourier
    @unpack valences = bulk_system.molsys.properties.monomers
    PsiC = fields.excess.PsiC[1]
    Rsys_star = CartesianIndices(f_hat)

    # Calculate excess chemical contribution by weighting point-wise derivative
    @. f_hat = surf_hat
    for j in axes(rho_beads_hat, ndims(rho_beads_hat))
        for K_star in Rsys_star
            f_hat[K_star] += rho_beads_hat[K_star, j] * valences[j]
        end

        # f_hat[1] += PsiC * valences[j]
    end

    @. f_hat *= weight_fft_K
    
    for j in axes(mu_ex_hat, ndims(mu_ex_hat))
        for K_star in Rsys_star
            mu_ex_hat[K_star, j] += f_hat[K_star] * valences[j]
        end
    end

    @unpack mu_ex_K = fields.excess
    for j in axes(mu_ex_K, ndims(mu_ex_K))
        for K in CartesianIndices(axes(mu_ex_K)[1:end-1])
            mu_ex_K[K, j] += PsiC * valences[j]
        end
    end
end