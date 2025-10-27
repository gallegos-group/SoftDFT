struct SquareWellFunctional{O} <: AbstractFunctional
    weight_fft_K :: Array{Float64, O}
    model        :: SquareWellFreeEnergy
end

"""
    construct_functional(::Type{SquareWellFunctional}, molsys::MolecularSystem, bulk::BulkState, geometry::CartesianCoord) -> SquareWellFunctional

Constructs a `SquareWellFunctional` for square-well interactions in DFT.

# Arguments
- `molsys`: Molecular system containing monomer diameters and interaction parameters.
- `bulk`: Bulk state containing the MeanFieldSW model.
- `geometry`: Spatial discretization of the system (`NP`, `bin_width`).

# Returns
- A fully initialized `SquareWellFunctional` object.
"""
function construct_functional(::Type{SquareWellFunctional}, molsys::MolecularSystem, bulk::BulkState, geometry::CartesianCoord)
    @unpack NP, bin_width, mirrored, offset = geometry

    @unpack monomers, pairs = molsys.properties
    @unpack diameters = monomers
    
    sequ_swpairs = pairs[:square_well][:pairs]
    @unpack lambdas = pairs[:square_well]

    NB = length(sequ_swpairs)

    # use calculation domain
    NP_star   = compute_full_domain(NP, mirrored, offset)
    dims_star = (NP_star..., NB)

    weight_fft_K = zeros(Float64, dims_star)
    weight_fft = zeros(Float64, NP_star)
    Rsys_star = CartesianIndices(NP_star)

    for (u, (i, j)) in enumerate(sequ_swpairs)
        λ = lambdas[u]

        # Outer well radius
        b_outer = λ * (diameters[i] + diameters[j]) / 2.0
        step_func_fft!(weight_fft, bin_width, b_outer)
        for K_star in Rsys_star
            weight_fft_K[K_star, u] = weight_fft[K_star]
        end

        # Inner well radius (excluded core)
        b_inner = (diameters[i] + diameters[j]) / 2.0
        step_func_fft!(weight_fft, bin_width, b_inner)
        for K_star in Rsys_star
            weight_fft_K[K_star, u] -= weight_fft[K_star]
        end
    end

    model = only(filter(e -> isa(e, SquareWellFreeEnergy), bulk.fe_model))
    return SquareWellFunctional{ndims(weight_fft_K)}(weight_fft_K, model)
end



function eval_fe_functional(bulk_system :: IsingLST, geometry :: CoordSystem, fields :: SpatialFields, fe_term :: SquareWellFunctional)
    @unpack mirrored, offset = geometry

    # Calculate excess fe contribution by weighting point-wise derivative
    @unpack species, pairs, energys = bulk_system.molsys.properties.pairs[:square_well]
    @unpack rho_beads_hat, plan_backward, f_hat = fields.fourier
    @unpack weight_fft_K = fe_term
    rho_beads_K = fields.rho_K.beads
    Rsys_star = CartesianIndices(f_hat)

    fe = 0.0
    for (_, i) in enumerate(species)
        @. f_hat = 0.0
        for (_, j) in enumerate(species)
            pair = i < j ? (i, j) : (j, i)
            index = findfirst(==(pair), pairs)

            if !isnothing(index)
                for K_star in Rsys_star
                    f_hat[K_star] += rho_beads_hat[K_star, j] * weight_fft_K[K_star, index] * energys[index]
                end
            end
        end

        f_hat = plan_backward * f_hat

        for K in CartesianIndices(axes(rho_beads_K)[1:end-1])
            K_star = to_star_index(K, offset, mirrored)
            fe -= 0.5 * rho_beads_K[K, i] * real(f_hat[K_star])
        end
    end

    return fe
end


function eval_mu_functional!(bulk_system :: IsingLST, fields :: SpatialFields, fe_term :: SquareWellFunctional)
    
    # Calculate excess chemical contribution by weighting point-wise derivative
    @unpack species, pairs, energys = bulk_system.molsys.properties.pairs[:square_well]
    @unpack rho_beads_hat, plan_backward, f_hat, mu_ex_hat = fields.fourier
    @unpack weight_fft_K = fe_term
    Rsys_star = CartesianIndices(f_hat)

    for (_, i) in enumerate(species)
        @. f_hat = 0.0
        for (_, j) in enumerate(species)
            pair = i < j ? (i, j) : (j, i)
            index = findfirst(==(pair), pairs)

            if !isnothing(index)
                for K_star in Rsys_star
                    f_hat[K_star] += rho_beads_hat[K_star, j] * weight_fft_K[K_star, index] * energys[index]
                end
            end
        end

        for K_star in Rsys_star
            mu_ex_hat[K_star, i] -= real(f_hat[K_star])
        end
    end
end