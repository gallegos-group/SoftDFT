function eval_spatial_density(bulk_system :: IsingLST, geometry :: CartesianCoord, fields :: SpatialFields)

    rho1_segments_K = [
            zeros(Float64, 
            size(fields.rho_K.segments[u])) 
            for u in eachindex(fields.rho_K.segments)]

    rho1_bonds_K = zeros(Float64, size(fields.rho_K.bonds))

    for (u, evaluator) in enumerate(bulk_system.bulk.evaluation)
        species_segments_K, species_bonds_K = 
            eval_spatial_density(u, bulk_system, geometry, fields, evaluator)
        
        for Kj in CartesianIndices(rho1_segments_K[u])
            rho1_segments_K[u][Kj] += species_segments_K[Kj]
        end

        for Kj in CartesianIndices(rho1_bonds_K)
            rho1_bonds_K[Kj] += species_bonds_K[Kj]
        end
    end
    
    # rho_simulated_K = fields.rho_K.simulated
    # for Kj in CartesianIndices(rho_simulated_K)
    #     rho1_beads_K[Kj] += rho_simulated_K[Kj]
    # end

    return rho1_segments_K, rho1_bonds_K
end

function eval_spatial_density(u, bulk_system, geometry, fields, evaluator :: SimulationEval)
    # Do nothing
end

function eval_spatial_density(u, bulk_system, geometry, fields, evaluator :: AnalyticalEval)
    @unpack species_model = bulk_system.molsys.configurations[u]

    return eval_spatial_density(u, bulk_system, geometry, fields, species_model)
end