function accumulate_densities(
        u,
        norm,
        gP,
        gC,
        bulk_system :: IsingLST, 
        geometry :: CartesianCoord, 
        fields :: SpatialFields, 
        species_model :: freely_jointed)

    # Related to species u
    config_u = bulk_system.molsys.configurations[u]
    mu_species_u = bulk_system.bulk.mu_species[u]

    species_bonds_K = zeros(Float64, size(fields.rho_K.bonds))
    species_segments_K = zeros(Float64, size(fields.rho_K.segments[u]))

    # Fields
    @unpack mu_ex_K, Ext, trapez = fields.excess

    # Monomers
    @unpack diameters, delta_muH = bulk_system.molsys.properties.monomers

    # Geometry
    @unpack bin_width, NP, features = geometry

    # Topology
    @unpack state_family, topology = config_u
    @unpack parents, children = topology
    @unpack bond_types = bulk_system.molsys.properties
    
    for parent in eachindex(state_family)
        for (idx_i, state_i) in enumerate(state_family[parent])
            for K in CartesianIndices(NP)
                contr = norm*exp(mu_species_u - mu_ex_K[K, state_i]  
                                    - Ext[K, state_i] - delta_muH[state_i]) *
                                    sum(@views gP[K, :, idx_i, parent])
            
                if isempty(children[parent])
                    temp = sum(@views gC[K, :, 1, idx_i, parent])
                else
                    temp = 0.0
                    for (idx_c, child) in enumerate(children[parent])

                        temp1 = 0.0
                        for (idx_j, state_j) in enumerate(state_family[child])
                            bond, _, _ = get_bond_info(state_i, state_j, bond_types)

                            temp2 = gC[K, idx_j, idx_c, idx_i, parent]
                            for (idx_oc, other_child) in enumerate(children[parent])
                                if other_child != child
                                    temp2 *= sum(@views gC[K, :, idx_oc, idx_i, parent])
                                end
                            end

                            species_bonds_K[K, 1, bond] += contr*temp2
                            temp1 += temp2
                        end
                        temp = temp1
                    end
                end

                species_segments_K[K, idx_i, parent] += contr*temp
            end
        end
    end

    for child in eachindex(state_family)
        if isempty(parents[child])
            continue  # root segment; no parent, so cannot contribute to bond direction 2
        end
        
        for (idx_i, state_i) in enumerate(state_family[child])
            for K in CartesianIndices(NP)
                contr = norm*exp(mu_species_u - mu_ex_K[K, state_i]  - Ext[K, state_i] - delta_muH[state_i])

                for (idx_gc, grandchildren) in enumerate(children[child])
                    contr *= sum(@views gC[K, :, idx_gc, idx_i, child])
                end

                for (_, parent) in enumerate(parents[child])
                    for (idx_j, state_j) in enumerate(state_family[parent])
                        bond, _, _ = get_bond_info(state_i, state_j, bond_types)

                        species_bonds_K[K, 2, bond] += contr*gP[K, idx_j, idx_i, child]
                    end
                end
            end
        end
    end

    return species_segments_K, species_bonds_K
end