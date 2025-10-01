function eval_fe_functional(bulk_system :: IsingLST, geometry :: CoordSystem, fields :: SpatialFields, functionals :: Vector{AbstractFunctional})

    fe = 0.0
    for (_, fe_term) in enumerate(functionals)
        fe += eval_fe_functional(bulk_system, geometry, fields, fe_term)
    end

    return fe
end