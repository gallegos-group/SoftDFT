
function process_fixed(dataset, molsys, fields, geometry)
    species_data = dataset["species"]
    species_properties = molsys.properties.species
    species_properties[:fixed_position] = [fill("", length(cfg.sequence)) for cfg in molsys.configurations]
    @unpack input_densities = species_properties
    
    for (u, fixed_struct) in enumerate(fields.fixed)
        if haskey(species_data[u], "fixed")
            fixed_data = species_data[u]["fixed"]
            sequence_len = length(molsys.configurations[u].sequence)

            # --- Check that bulk density is zero ---
            if input_densities[u] > 0.0
                error("Species $u defines a fixed segment but also has nonzero input density ($(input_densities[u])).\n" *
                      "Simultaneous bulk and brush densities are not currently supported.")
            end

            # --- Required: fixed density ---
            if haskey(fixed_data, "density")
                fixed_struct.density[1] = fixed_data["density"]
            else
                error("Species $u contains 'fixed' but does not specify a density.")
            end

            # --- Required: fixed segments and positions ---
            if haskey(fixed_data, "segment") && haskey(fixed_data, "position")
                segment_data = fixed_data["segment"]
                position_data = fixed_data["position"]

                if length(segment_data) != length(position_data)
                    error("Mismatch in number of fixed segments and specified positions for species $u.")
                end

                for (i, L) in enumerate(segment_data)
                    fixed_struct.segments[L] = true
                    species_properties[:fixed_position][u][L] = position_data[i]
                end
            else
                error("Species $u contains 'fixed' but is missing 'segment' or 'position'.")
            end
        end
    end

    determine_fixed(molsys, fields, geometry)
end