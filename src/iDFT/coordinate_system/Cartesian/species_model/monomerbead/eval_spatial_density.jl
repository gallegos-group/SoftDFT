function eval_spatial_density(
        u, 
        bulk_system, 
        geometry :: CartesianCoord, 
        fields :: SpatialFields, 
        species_model :: monomerbead)

        norm = handle_fixed(u, bulk_system, geometry, fields, species_model)
        return accumulate_densities(u, norm, bulk_system, geometry, fields, species_model)
end