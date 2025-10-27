function determine_fixed(fixed_position, molsys, fields, geometry::CartesianZ)

    @unpack diameters = molsys.properties.monomers
    @unpack bin_width, offset, dimensions = geometry

    idx_LW, idx_RW = compute_wall_indices(dimensions, bin_width)

    for (u, fixed_struct) in enumerate(fields.fixed)
        config = molsys.configurations[u]
        sequence = config.sequence
        pos_tags = fixed_position[u]
        coords = fixed_struct.coordinates

        for (L, isfixed) in enumerate(fixed_struct.segments)
            if isfixed
                j = sequence[L]
                dz = round(Int, diameters[j] / bin_width[end] / 2.0)

                if pos_tags[L] == "left"
                    coords[L] = (idx_LW[end] + dz,)
                elseif pos_tags[L] == "right"
                    coords[L] = (idx_RW[end] - dz,)
                else
                    error("Invalid grafting position tag for species $u, segment $L: '$(pos_tags[L])'")
                end
            end
        end
    end
end

function determine_fixed(fixed_position, molsys, fields, geometry::CartesianXYZ)

    @unpack diameters = molsys.properties.monomers
    @unpack NP, bin_width, offset, dimensions = geometry

    idx_LW, idx_RW = compute_wall_indices(dimensions, bin_width)

    z_count = 0
    left_count = 0
    right_count = 0

    for (u, fixed_struct) in enumerate(fields.fixed)
        config = molsys.configurations[u]
        sequence = config.sequence
        pos_tags = fixed_position[u]
        coords = fixed_struct.coordinates

        if count(fixed_struct.segments) > 1
            error("Only one fixed segment per chain is allowed in 3D (species $u).")
        end

        for (L, isfixed) in enumerate(fixed_struct.segments)
            if isfixed
                j = sequence[L]
                dz = round(Int, diameters[j] / bin_width[end] / 2.0)
                x_center = round(Int, NP[1] / 2) + 1
                y_center = round(Int, NP[2] / 2) + 1

                if dimensions[1] != dimensions[2]
                    error("For 3D brushes, wall dimensions in x and y must match.")
                end

                # Normalize fixed density based on wall area
                println("Grafting density is currently specified as $(fixed_struct.density[1])")
                fixed_struct.density[1] = 1.0 / (dimensions[1] * dimensions[2])
                println("Updated grafting density: $(fixed_struct.density[1]) based off system dimensions: x = $(dimensions[1]) and y = $(dimensions[2]).")

                if pos_tags[L] == "left"
                    right_count += 1
                    coords[L] = (x_center, y_center, idx_LW[end] + dz)
                elseif pos_tags[L] == "right"
                    left_count += 1
                    coords[L] = (x_center, y_center, idx_RW[end] - dz)
                else
                    error("Invalid grafting position tag for species $u, segment $L: '$(pos_tags[L])'")
                end

                z_count += 1
            end
        end
    end

    if left_count > 1 || right_count > 1
        error("Only one brush may be attached to each wall in 3D.")
    end
end

function compute_wall_indices(dimensions, bin_width)
    idx_LW = round.(Int, 0.0 ./ bin_width) .+ 1
    idx_RW = round.(Int, dimensions ./ bin_width) .+ 1
    return idx_LW, idx_RW
end
