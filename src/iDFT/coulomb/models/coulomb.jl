"""
    struct coulFunctional{N} <: AbstractFunctional

Inhomogeneous DFT functional for the Coulomb contribution.

# Type Parameters
- `N`: Dimensionality of the simulation grid.

# Fields
- `weight_fft_K::Array{Float64,N}`: Fourier-space Coulomb kernel.
- `model::CoulombFreeEnergy`: Precomputed bulk free energy model.
"""
struct CoulombFFT{N} <: AbstractCoulomb
    weight_fft_K :: Array{Float64, N}
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
function construct_coulomb( :: Type{CoulombFFT},
                              molsys::MolecularSystem,
                              bulk::BulkState,
                              geometry::CartesianCoord)
    @unpack NP, bin_width, mirrored, offset = geometry
    @unpack bjerrum_length = molsys.properties.system

    # work on the calculation domain
    NP_star = compute_full_domain(NP, mirrored, offset)

    weight_fft_K = zeros(Float64, NP_star)
    coul_func_fft!(weight_fft_K, bin_width, bjerrum_length)

    return CoulombFFT{ndims(weight_fft_K)}(weight_fft_K)
end




function eval_coulomb_energy(bulk_system :: IsingLST, geometry :: CoordSystem, fields :: SpatialFields, coulomb :: CoulombFFT)
    @unpack mirrored, offset = geometry

    @unpack weight_fft_K = coulomb
    @unpack surf_hat, rho_beads_hat, f_hat, mu_ex_hat, plan_backward = fields.fourier
    @unpack valences = bulk_system.molsys.properties.monomers
    PsiC = fields.excess.PsiC[1]
    Psi = fields.excess.Psi

    Rsys_star = CartesianIndices(f_hat)

    @. f_hat = surf_hat
    for j in axes(rho_beads_hat, ndims(rho_beads_hat))
        for K_star in Rsys_star
            f_hat[K_star] += rho_beads_hat[K_star, j] * valences[j]
        end
    end
    
    f_hat = plan_backward * f_hat

    fe = 0.0
    for K in CartesianIndices(Psi)
        K_star = to_star_index(K, offset, mirrored)
        fe += real(f_hat[K_star]) * (Psi[K] + PsiC)/ 2.0
    end

    return fe
end


"""
    solve_coulomb!(dft_system)

Updates `fields.excess.Psi` and `fields.excess.PsiC`.
"""
function solve_coulomb(bulk_system, geometry, fields, coulomb::CoulombFFT)

    # Step 1: Solve Poisson
    Psi = solve_poisson_eq(bulk_system, geometry, fields, coulomb)

    # Step 2: Enforce Charge Neutrality
    PsiC = compute_neutrality_shift(bulk_system, geometry, fields)

    return Psi, PsiC
end


function solve_poisson_eq(bulk_system, geometry, fields, coulomb :: CoulombFFT)
    @unpack mirrored, offset = geometry

    @unpack weight_fft_K = coulomb
    @unpack surf_hat, rho_beads_hat, f_hat, plan_backward = fields.fourier
    @unpack valences = bulk_system.molsys.properties.monomers

    Psi = zeros(Float64, size(fields.excess.Psi))

    Rsys_star = CartesianIndices(f_hat)

    # Calculate excess chemical contribution by weighting point-wise derivative
    @. f_hat = surf_hat
    for j in axes(rho_beads_hat, ndims(rho_beads_hat))
        for K_star in Rsys_star
            f_hat[K_star] += rho_beads_hat[K_star, j] * valences[j]
        end
    end

    @. f_hat *= weight_fft_K
    
    f_hat = plan_backward * f_hat

    for K in CartesianIndices(Psi)
        K_star = to_star_index(K, offset, mirrored)
        Psi[K] = real(f_hat[K_star])
    end
    
    return Psi
end