function brush_height(dft_system; filename::String = "brush_height.txt")
    rho_segments   = dft_system.fields.rho_K.segments
    trapez         = dft_system.fields.excess.trapez
    molsys         = dft_system.bulk_system.molsys
    fixed_species  = dft_system.fields.fixed
    species_names  = molsys.properties.species[:species]
    species_configs = molsys.configurations
    @unpack NP, bin_width = dft_system.geometry

    Rsys = CartesianIndices(NP)

    # --- NEW: skip output entirely if no species has any fixed segments
    has_fixed_any = any(i -> any(fixed_species[i].segments), eachindex(species_names))
    if !has_fixed_any
        return nothing
    end

    # Precompute z along the normal direction to avoid recomputing in inner loops
    z_vec = [(k - 1) * bin_width[end] for k in 1:NP[end]]

    open(filename, "w") do io
        println(io, "# Brush heights per species (computed as 2 * <z>)")
        println(io, "# Format: species brush_height")

        for (i, species_name) in enumerate(species_names)
            config = species_configs[i]
            if !any(fixed_species[i].segments)
                continue  # only report species with at least one fixed segment
            end

            rho_segments_i = rho_segments[i]
            height_numer = 0.0
            height_denom = 0.0
            n_segments = length(config.sequence)

            for seg in 1:n_segments
                for (idx, state_i) in enumerate(config.state_family[seg])
                    if fixed_species[i].segments[seg]
                        # Fixed segments: use ρ/bin_width (line density per z-bin)
                        for K in Rsys
                            k1 = K[1]
                            val = rho_segments_i[K, idx, seg] / bin_width[end]
                            height_numer += val * z_vec[k1]
                            height_denom += val
                        end
                    else
                        # Free segments: use trapezoidal weights
                        for K in Rsys
                            k1 = K[1]
                            val = rho_segments_i[K, idx, seg] * trapez[K, state_i]
                            height_numer += val * z_vec[k1]
                            height_denom += val
                        end
                    end
                end
            end

            height = height_denom > 0 ? 2 * height_numer / height_denom : 0.0
            println(io, "$(species_name) $(height)")
        end
    end

    return nothing
end