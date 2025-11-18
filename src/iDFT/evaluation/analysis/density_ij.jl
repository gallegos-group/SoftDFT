
function calc_density_ij(u, segment_i, coord_i, dft_system)

    @unpack bulk_system, geometry, fields = dft_system

    @unpack species_model = bulk_system.molsys.configurations[u]
    
    # Geometry
    @unpack NP, bin_width, features = geometry

    # Fixed Segments
    fixed_species_u = fields.fixed[u]
    fixed_segments = fixed_species_u.segments
    fixed_coordinates = fixed_species_u.coordinates[1]
    fields.fixed[u].density[1] = 1.0

    fixed_segments[segment_i] = true
    fixed_coordinates[segment_i] = (floor(Int, coord_i/bin_width[1]) + 1,) 

    species_segments_K, _ = eval_spatial_density(u, bulk_system, geometry, fields, species_model)

    write_field_array(species_segments_K[:, 1, :], "density_ij.txt", geometry.NP)
    
    fields.fixed[u].density[1] = 0.0
    fixed_segments[segment_i] = false
    fixed_coordinates[segment_i] = () 
end