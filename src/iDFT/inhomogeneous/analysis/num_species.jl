function get_num_species(
        bulk_system :: IsingLST, 
        geometry :: CartesianCoord, 
        fields :: SpatialFields)

    @unpack configurations = bulk_system.molsys
    @unpack valences = bulk_system.molsys.properties.monomers
    @unpack NP, bin_width, features = geometry
    @unpack fixed = fields
    @unpack trapez = fields.excess

    rho_segments_K = fields.rho_K.segments
    bin = prod(bin_width)

    Rsys = CartesianIndices(NP)

    Num_species = zeros(Float64, length(configurations))
    for (u, config) in enumerate(configurations)
        rho_segments_u = rho_segments_K[u]
        for seg_i in eachindex(config.sequence)
            for (idx_i, state_i) in enumerate(config.state_family[seg_i])
                if fixed[u].segments[seg_i]
                    K = CartesianIndex(fixed[u].coordinates[1][seg_i])
                    num = rho_segments_u[K, idx_i, seg_i]
                else
                    num = 0.0
                    for K in Rsys
                        num += rho_segments_u[K, idx_i, seg_i]*trapez[K, state_i]*bin
                    end
                end

                Num_species[u] += num / length(config.sequence)
            end
        end
    end

    return Num_species
end