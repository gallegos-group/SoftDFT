function accumulate_densities(
        u,
        norm,
        bulk_system :: IsingLST, 
        geometry :: CartesianCoord, 
        fields :: SpatialFields, 
        species_model :: monomerbead)

    # Related to species u
    config_u = bulk_system.molsys.configurations[u]
    mu_species_u = bulk_system.bulk.mu_species[u]
    species_segments_K = zeros(Float64, size(fields.rho_K.segments[u]))

    # Fields
    @unpack mu_ex_K, Ext, trapez = fields.excess

    # Monomers
    @unpack diameters, delta_muH = bulk_system.molsys.properties.monomers

    # Geometry
    @unpack bin_width, NP = geometry

    # Topology
    @unpack state_family, topology = config_u
    @unpack parents, children = topology
    @unpack bond_types = bulk_system.molsys.properties

    # Fixed Segments
    fixed_segments = fields.fixed[u].segments

    for parent in eachindex(state_family)
        for (idx_i, state_i) in enumerate(state_family[parent])
            for K in CartesianIndices(NP)
                contr = norm*exp(mu_species_u - mu_ex_K[K, state_i]  
                                    - Ext[K, state_i] - delta_muH[state_i])
    
                species_segments_K[K, idx_i, parent] += contr
            end
        end
    end

    return species_segments_K, zeros(Float64, size(fields.rho_K.bonds))
end