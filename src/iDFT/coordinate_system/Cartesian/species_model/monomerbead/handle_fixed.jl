function handle_fixed(
        u,
        bulk_system :: IsingLST, 
        geometry :: CartesianCoord, 
        fields :: SpatialFields, 
        species_model :: monomerbead)

    # Fields
    @unpack mu_ex_K, Ext, trapez = fields.excess

    # Monomers
    @unpack delta_muH = bulk_system.molsys.properties.monomers

    # Geometry
    @unpack bin_width, NP, features = geometry

    # Topology
    @unpack state_family = bulk_system.molsys.configurations[u]

    # Fixed Segments
    fixed_density = fields.fixed[u].density[1]
    fixed_segments = fields.fixed[u].segments
    fixed_coordinates = fields.fixed[u].coordinates

    if any(fixed_segments)
        norm = fixed_density
        bulk_system.bulk.mu_species[u] = 0.0
        lambda = 0.0

        parent = findfirst(fixed_segments)
        @assert parent !== nothing "No fixed segment available for lambda normalization."

        K = fixed_coordinates[1][parent]
            for (idx_i, state_i) in enumerate(state_family[parent])

                contr = exp(-mu_ex_K[K, state_i] - Ext[K, state_i] - delta_muH[state_i])
                lambda += contr
            end
    else
        lambda = 1.0
        norm = 1.0
    end
    fields.fixed[u].lambda[1] = lambda

    if !isfinite(lambda) || lambda <= 0
        error("Normalization error: λ = $lambda for species $u. Check propagator integrity.")
    end
    norm /= lambda

    return norm
end