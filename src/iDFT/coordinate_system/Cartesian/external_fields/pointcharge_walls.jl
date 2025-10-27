#   external_field:
    # pointcharged_walls:
    #   dims: [1]
    #   charge: [0.0, -0.0]

function eval_External(molsys, geometry::CartesianCoord, fields, ext_field :: ExtPointChargeWalls)
    @unpack NP, bin_width, offset, mirrored = geometry

    n_dims = length(NP)
    inv_voxel_volume = 1.0 / prod(bin_width)

    @unpack surf_hat = fields.fourier

    charges  = ext_field.charges
    walls    = ext_field.positions

    Rsys = CartesianIndices(NP)  # true domain

    for K in Rsys
        idx = Tuple(K)

        for (v, wall_positions) in enumerate(walls)
            isempty(wall_positions) && continue

            # Check if idx lies in the wall plane (fix all  except v)
            in_plane = true
            for v1 in 1:n_dims
                if v1 != v
                    center_idx = round(Int, NP[v1] / 2) + 1
                    if idx[v1] != center_idx
                        in_plane = false
                        break
                    end
                end
            end
            in_plane || continue

            @inbounds for (i, wall_pos) in enumerate(wall_positions)
                grid_idx = round(Int, wall_pos / bin_width[v]) + 1
                K_star = to_star_index(K, offset, mirrored)
                if idx[v] == grid_idx
                    surf_hat[K_star] += charges[v][i] * inv_voxel_volume
                end
            end
        end
    end

    apply_mirroring!(surf_hat, geometry)
end