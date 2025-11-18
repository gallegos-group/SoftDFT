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

@generated function _eval_spatial_density_tuple!(rho1_segments_K, rho1_bonds_K,
                                                 bulk_system::IsingLST,
                                                 geometry::CartesianCoord,
                                                 fields::SpatialFields,
                                                 evals::Tuple{Vararg{AbstractEvaluation, N}}) where {N}
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

    
    get_onebody_potential!(bulk_system, fields)
    
    _eval_spatial_density_tuple!(rho1_segments_K, rho1_bonds_K,
                                 bulk_system, geometry, fields,
                                 bulk_system.bulk.evaluation)
    return rho1_segments_K, rho1_bonds_K
end





# in place

function eval_spatial_density!(u, bulk_system, geometry, fields, evaluator :: SimulationEval)
end

function eval_spatial_density!(u, bulk_system, geometry, fields, evaluator :: AnalyticalEval)
    @unpack species_model = bulk_system.molsys.configurations[u]

    return eval_spatial_density!(u, bulk_system, geometry, fields, species_model)
end


# === Recursive helper ===
@generated function _eval_spatial_density_tuple!(bulk_system::IsingLST,
                                                 geometry::CartesianCoord,
                                                 fields::SpatialFields,
                                                 evals::NTuple{N, AnalyticalEval}) where {N}
    quote
        Base.Cartesian.@nexprs $N i -> begin
            eval_spatial_density!(i, bulk_system, geometry, fields, evals[i])
        end
    end
end

@generated function _eval_spatial_density_tuple!(bulk_system::IsingLST,
                                                 geometry::CartesianCoord,
                                                 fields::SpatialFields,
                                                 evals::Tuple{Vararg{AbstractEvaluation, N}}) where {N}
    quote
        Base.Cartesian.@nexprs $N i -> begin
            eval_spatial_density!(i, bulk_system, geometry, fields, evals[i])
        end
    end
end


function eval_spatial_density!(bulk_system::IsingLST,
                              geometry::CartesianCoord,
                              fields::SpatialFields)

    for (u, segments) in enumerate(fields.rho_K.segments)
        @. segments = 0.0
    end
    @. fields.rho_K.bonds = 0.0

    get_onebody_potential!(bulk_system, fields)

    _eval_spatial_density_tuple!(bulk_system, geometry, fields, bulk_system.bulk.evaluation)

end

function get_onebody_potential!(bulk_system, fields)

    @unpack mu_ex_K, Ext, Psi, PsiC, lambda_K = fields.excess
    @unpack delta_muH, valences = bulk_system.molsys.properties.monomers

    for K in CartesianIndices(Psi)
        ψ = Psi[K] + PsiC[1]
        for j in eachindex(valences)
            lambda_K[K, j] = mu_ex_K[K, j] + ψ * valences[j] + Ext[K, j] + delta_muH[j]
        end
    end
end