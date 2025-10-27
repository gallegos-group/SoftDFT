
#   external_field:
#     hard_walls:
#       dims: [1]

function eval_External(molsys, geometry::CartesianCoord, fields, ext_field :: ExtHardWalls)
    @unpack Ext, trapez = fields.excess
    @unpack NP, bin_width, offset, mirrored = geometry
    @unpack diameters = molsys.properties.monomers

    Rsys = CartesianIndices(NP)

    walls = ext_field.positions

    for K in Rsys
        idx = Tuple(K)

        for (v, wall_positions) in enumerate(walls)
            isempty(wall_positions) && continue

            wall_lo = wall_positions[1]
            wall_hi = wall_positions[2]

            for (j, d) in enumerate(diameters)
                lo_idx = round(Int, (wall_lo + d / 2.0) / bin_width[v]) + 1
                hi_idx = round(Int, (wall_hi - d / 2.0) / bin_width[v]) + 1

                if idx[v] < lo_idx
                    Ext[K, j] += 1.0e20
                    trapez[K, j] = 1.0
                elseif idx[v] == lo_idx
                    Ext[K, j] += 0.0
                    trapez[K, j] = 0.5
                end

                if idx[v] > hi_idx
                    Ext[K, j] += 1.0e20
                    trapez[K, j] = 1.0
                elseif idx[v] == hi_idx
                    Ext[K, j] += 0.0
                    trapez[K, j] = 0.5
                end
            end
        end
    end
end


function contact_value_theorem(molsys, geometry::CartesianCoord, fields, ext_field :: ExtHardWalls)

    @unpack Ext, trapez = fields.excess
    @unpack NP, bin_width, offset, mirrored = geometry
    @unpack diameters = molsys.properties.monomers

    walls = ext_field.positions

    rho_beads_K = fields.rho_K.beads

    contr = zeros(Float64, (2,length(geometry.NP)))
    for (v, wall_positions) in enumerate(walls)
        isempty(wall_positions) && continue

        wall_lo = wall_positions[1]
        wall_hi = wall_positions[2]

        for (j, d) in enumerate(diameters)
            lo_idx = round(Int, (wall_lo + d / 2.0) / bin_width[v]) + 1
            hi_idx = round(Int, (wall_hi - d / 2.0) / bin_width[v]) + 1

            K = CartesianIndex((NP[1:end-1]..., lo_idx))
            contr[1, v] += rho_beads_K[K, j]

            K = CartesianIndex((NP[1:end-1]..., hi_idx))
            contr[2, v] += rho_beads_K[K, j]               
        end
    end

    return contr
end