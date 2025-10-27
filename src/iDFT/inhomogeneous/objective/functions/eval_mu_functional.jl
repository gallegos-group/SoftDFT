function eval_mu_functionals!(
    bulk_system :: IsingLST,
    geometry    :: CoordSystem,
    fields      :: SpatialFields,
    functionals :: F
        ) where {F <: Tuple{Vararg{AbstractFunctional}}}

    @unpack NP, mirrored, offset = geometry

    @unpack mu_ex_K, lng_K = fields.excess
    @unpack mu_ex_hat, mu_sim_hat, lng_hat, f_hat, plan_backward = fields.fourier

    @. mu_ex_K = 0.0
    @. lng_K = 0.0

    @. mu_ex_hat = mu_sim_hat  # Zero if no external field simulation
    @. lng_hat = 0.0

    eval_mu_functional!(bulk_system, fields, functionals)

    Rsys_star = CartesianIndices(f_hat)
    Rsys = CartesianIndices(NP)

    for j in axes(mu_ex_K, ndims(mu_ex_K))
        for K_star in Rsys_star
            f_hat[K_star] = mu_ex_hat[K_star, j]
        end

        f_hat = plan_backward * f_hat

        for K in Rsys
            K_star = to_star_index(K, offset, mirrored)
            mu_ex_K[K, j] += real(f_hat[K_star])
        end
    end

    for u in axes(lng_K, ndims(lng_K))
        for idx in axes(lng_K, ndims(lng_K)-1)
            for K_star in Rsys_star
                f_hat[K_star] = lng_hat[K_star, idx, u]
            end

            f_hat = plan_backward * f_hat

            for K in Rsys
                K_star = to_star_index(K, offset, mirrored)
                lng_K[K, idx, u] += real(f_hat[K_star])
            end
        end
    end
end


# Base case: empty tuple → stop
@inline function eval_mu_functional!(
    bulk_system::IsingLST,
    fields::SpatialFields,
    ::Tuple{}
)
    return
end

# Recursive case: process the first functional, then recurse on the rest
@inline function eval_mu_functional!(
    bulk_system::IsingLST,
    fields::SpatialFields,
    functionals::Tuple
)
    eval_mu_functional!(bulk_system, fields, functionals[1])
    eval_mu_functional!(bulk_system, fields, Base.tail(functionals))
end