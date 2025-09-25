function process_guess(dft_system)
    return process_guess(dft_system.bulk_system, dft_system.geometry, dft_system.fields)
end

function process_guess(bulk_system, geometry, fields)
    @unpack bulk, molsys = bulk_system

    rho_segments_K = fields.rho_K.segments
    rho_beads_K = fields.rho_K.beads
    rho_bonds_K = fields.rho_K.bonds

    @unpack Ext, mu_ex_K, lng_K = fields.excess

    rho_segments = bulk.rho.segments
    rho_beads = bulk.rho.beads
    rho_bonds = bulk.rho.bonds
    sequ_bonds = molsys.properties.bond_types

    Rsys = CartesianIndices(geometry.NP)

    for (u, config) in enumerate(molsys.configurations)
        for seg_i in eachindex(config.sequence)                
            for (idx_i, state_i) in enumerate(config.state_family[seg_i])
                for K in Rsys
                    if Ext[K, state_i] < 1e2
                        rho_segments_K[u][K, idx_i, seg_i] = 
                                rho_segments[u][idx_i, seg_i] * exp(-Ext[K, state_i])
                    end
                end
            end
        end
    end

    for j in axes(rho_beads_K, ndims(rho_beads_K))
        for K in Rsys
            mu_ex_K[K, j] = bulk.mu_ex[j]
            if Ext[K, j] < 1e2
                rho_beads_K[K, j] = rho_beads[j]
            end
        end
    end
    
    for u in axes(rho_bonds_K, ndims(rho_bonds_K))
        for idx in axes(rho_bonds_K, ndims(rho_bonds_K)-1)
            for K in Rsys
                lng_K[K, idx, u] = bulk.lng[u]

                if Ext[K, sequ_bonds[u][idx]] < 1e2
                    rho_bonds_K[K, idx, u] = rho_bonds[idx, u]
                end
            end
        end
    end
end