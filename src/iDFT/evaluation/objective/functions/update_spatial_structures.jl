function update_spatial_structures!(
    bulk_system :: IsingLST,
    geometry :: CartesianCoord, 
    fields :: SpatialFields)

    @unpack configurations = bulk_system.molsys
    @unpack fixed = fields
    @unpack NP, bin_width = geometry
    @unpack trapez = fields.excess
    
    rho_segments_K = fields.rho_K.segments
    rho_beads_K = fields.rho_K.beads
    
    bin = prod(bin_width)

    Rsys = CartesianIndices(NP)

    @. rho_beads_K = 0.0
    for (u, config) in enumerate(configurations)
        for seg_i in eachindex(config.sequence)                
            for (idx_i, state_i) in enumerate(config.state_family[seg_i])
                
                if fixed[u].segments[seg_i]
                    K = CartesianIndex(fixed[u].coordinates[1][seg_i])
                    rho_beads_K[K, state_i] += rho_segments_K[u][K, idx_i, seg_i]/bin/trapez[K, state_i]
                else
                    for K in Rsys
                        rho_beads_K[K, state_i] += rho_segments_K[u][K, idx_i, seg_i]
                    end
                end
            end
        end
    end
end