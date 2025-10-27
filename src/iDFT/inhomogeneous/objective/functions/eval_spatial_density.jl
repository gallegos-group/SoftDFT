# function eval_spatial_density(bulk_system :: IsingLST, geometry :: CartesianCoord, fields :: SpatialFields)

#     rho1_segments_K = [
#             zeros(Float64, 
#             size(fields.rho_K.segments[u])) 
#             for u in eachindex(fields.rho_K.segments)]

#     rho1_bonds_K = zeros(Float64, size(fields.rho_K.bonds))

#     for (u, evaluator) in enumerate(bulk_system.bulk.evaluation)
#         species_segments_K, species_bonds_K = 
#             eval_spatial_density(u, bulk_system, geometry, fields, evaluator)::Tuple{typeof(rho1_segments_K[u]), typeof(rho1_bonds_K)}
        
#         for Kj in CartesianIndices(rho1_segments_K[u])
#             rho1_segments_K[u][Kj] += species_segments_K[Kj]
#         end

#         for Kj in CartesianIndices(rho1_bonds_K)
#             rho1_bonds_K[Kj] += species_bonds_K[Kj]
#         end
#     end

#     return rho1_segments_K, rho1_bonds_K
# end

function eval_spatial_density(u, bulk_system, geometry, fields, evaluator :: SimulationEval)
    return copy(fields.rho_K.segments[u]), zeros(Float64, size(fields.rho_K.bonds))
end

function eval_spatial_density(u, bulk_system, geometry, fields, evaluator :: AnalyticalEval)
    @unpack species_model = bulk_system.molsys.configurations[u]

    return eval_spatial_density(u, bulk_system, geometry, fields, species_model)
end


# === Recursive helper ===
@generated function _eval_spatial_density_tuple!(rho1_segments_K, rho1_bonds_K,
                                                 bulk_system::IsingLST,
                                                 geometry::CartesianCoord,
                                                 fields::SpatialFields,
                                                 evals::NTuple{N, AnalyticalEval}) where {N}
    quote
        Base.Cartesian.@nexprs $N i -> begin
            segs, bonds = eval_spatial_density(i, bulk_system, geometry, fields, evals[i])
            @. rho1_segments_K[i] += segs
            @. rho1_bonds_K += bonds
        end
    end
end

function eval_spatial_density(bulk_system::IsingLST,
                              geometry::CartesianCoord,
                              fields::SpatialFields)
    rho1_segments_K = [zeros(Float64, size(fields.rho_K.segments[u]))
                       for u in eachindex(fields.rho_K.segments)]
    rho1_bonds_K = zeros(Float64, size(fields.rho_K.bonds))
    _eval_spatial_density_tuple!(rho1_segments_K, rho1_bonds_K,
                                 bulk_system, geometry, fields,
                                 bulk_system.bulk.evaluation)
    return rho1_segments_K, rho1_bonds_K
end