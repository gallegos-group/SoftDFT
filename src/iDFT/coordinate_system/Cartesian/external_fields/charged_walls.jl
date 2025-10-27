    
#   external_field:
    # charged_walls:
    #   dims: [1]
    #   charge: [0.0, -0.0]

function eval_External(molsys, geometry::CartesianCoord, fields, ext_field :: ExtChargedWalls)
    @unpack NP, bin_width, mirrored, offset = geometry

    n_dims     = length(NP)
    Rsys = CartesianIndices(NP)

    @unpack surf_hat = fields.fourier

    charges  = ext_field.charges
    walls    = ext_field.positions

    for K in Rsys
        idx = Tuple(K) # True domain indices

        for v in 1:n_dims
            isempty(walls[v]) && continue

            # Check if the index is in the plane orthogonal to v
            in_wall_plane = true
            invbin_area = 1.0
            for v1 in 1:n_dims
                if v1 == v
                    invbin_area /= bin_width[v1]
                else
                    invbin_area *= bin_width[v1]

                    lo = 1
                    hi = NP[v1]

                    if idx[v1] < lo || idx[v1] > hi
                        in_wall_plane = false
                        break
                    end
                end
            end

            if in_wall_plane
                for (i, wall_pos) in enumerate(walls[v])
                    grid_pos = round(Int, wall_pos / bin_width[v]) + 1
                    K_star = to_star_index(K, offset, mirrored)
                    if idx[v] == grid_pos
                        surf_hat[K_star] += charges[v][i] * invbin_area
                    end
                end
            end
        end
    end
    
    # Fill mirrored halves
    apply_mirroring!(surf_hat, geometry)
end


function contact_value_theorem(molsys, geometry::CartesianCoord, fields, ext_field :: ExtChargedWalls)
    @unpack NP, bin_width, mirrored, offset = geometry

    n_dims     = length(NP)

    @unpack surf_hat = fields.fourier

    charges  = ext_field.charges
    walls    = ext_field.positions

    @unpack bjerrum_length = molsys.properties.system

    contr = zeros(Float64, (2,length(geometry.NP)))
    for v in 1:n_dims
        isempty(walls[v]) && continue

        contr[1, v] = -2.0*pi*bjerrum_length*charges[v][1]^2
        contr[2, v] = -2.0*pi*bjerrum_length*charges[v][2]^2
    end
    
    return contr
end