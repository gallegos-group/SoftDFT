function eval_spatial_density(
        u, 
        bulk_system, 
        geometry :: CartesianCoord, 
        fields :: SpatialFields, 
        species_model :: freely_jointed)

        gC = compute_child_propagators(u, bulk_system, geometry, fields, species_model)
        gP = compute_parent_propagators(u, gC, bulk_system, geometry, fields, species_model)
        norm = handle_fixed(u, gP, gC, bulk_system, geometry, fields, species_model)
        return accumulate_densities(u, norm, gP, gC, bulk_system, geometry, fields, species_model)
end

function eval_spatial_density!(
        u, 
        bulk_system, 
        geometry :: CartesianCoord, 
        fields :: SpatialFields, 
        species_model :: freely_jointed)

        gC = compute_child_propagators(u, bulk_system, geometry, fields, species_model)
        gP = compute_parent_propagators(u, gC, bulk_system, geometry, fields, species_model)
        norm = handle_fixed(u, gP, gC, bulk_system, geometry, fields, species_model)
        return accumulate_densities!(u, norm, gP, gC, bulk_system, geometry, fields, species_model)
end