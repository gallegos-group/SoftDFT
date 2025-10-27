#   external_field:
    # sw_walls:
    #   dims: [1]
    #   monomers:
    #     B:
    #       energy: [[-3.0, -0.5]]
    #       lambda: [[1.5, 1.5]]

        # C:
        #   energy: [[-0.0, 0.0]]
        #   lambda: [[1.5, 1.5]]

function eval_External(molsys, geometry::CartesianCoord, fields, ext_field :: ExtSquareWellWalls)

    @unpack Ext, trapez = fields.excess
    @unpack NP, bin_width, offset, mirrored = geometry
    @unpack diameters = molsys.properties.monomers
    
    energys = ext_field.energys
    lambdas = ext_field.lambdas
    wall_positions = ext_field.positions
    ext_names = ext_field.names

    # Build lookup: external field name → index
    ext_name_to_idx = Dict(name => i for (i, name) in enumerate(ext_names))

    Rsys = CartesianIndices(NP)

    for K in Rsys
        idx = Tuple(K)

        for (dim, walls) in enumerate(wall_positions)
            isempty(walls) && continue

            wall_lo = walls[1]
            wall_hi = walls[2]

            for (mon_idx, mon_name) in enumerate(molsys.properties.species.monomers)
                # --- Skip if this monomer has no external field defined ---
                ext_idx = get(ext_name_to_idx, mon_name, nothing)
                ext_idx === nothing && continue

                d = diameters[mon_idx]
                λ = lambdas[dim][ext_idx]
                ε = energys[dim][ext_idx]

                lo_idx = round(Int, (wall_lo + λ[1] * d / 2) / bin_width[dim]) + 1
                hi_idx = round(Int, (wall_hi - λ[2] * d / 2) / bin_width[dim]) + 1

                if idx[dim] <= lo_idx
                    Ext[K, mon_idx] += ε[1]
                elseif idx[dim] >= hi_idx
                    Ext[K, mon_idx] += ε[2]
                end
            end
        end
    end
end


function contact_value_theorem(molsys, geometry::CartesianCoord, fields, ext_field :: ExtSquareWellWalls)
    @unpack Ext = fields.excess
    @unpack NP, bin_width, offset, mirrored = geometry
    @unpack diameters = molsys.properties.monomers

    # External field data
    energys = ext_field.energys
    lambdas = ext_field.lambdas
    wall_positions = ext_field.positions
    ext_names = ext_field.names

    # --- Map from external-field name → index ---
    ext_name_to_idx = Dict(name => i for (i, name) in enumerate(ext_names))

    rho_beads_K = fields.rho_K.beads
    
    contr = zeros(Float64, (2,length(geometry.NP)))
    for (dim, walls) in enumerate(wall_positions)
        isempty(walls) && continue

        wall_lo = walls[1]
        wall_hi = walls[2]

        for (mon_idx, mon_name) in enumerate(molsys.properties.species.monomers)
            ext_idx = get(ext_name_to_idx, mon_name, nothing)
            ext_idx === nothing && continue  # skip monomers not in external field

            d = diameters[mon_idx]
            λ = lambdas[dim][ext_idx]
            ε = energys[dim][ext_idx]

            lo_idx = round(Int, (wall_lo + λ[1] * d / 2) / bin_width[dim]) + 1
            hi_idx = round(Int, (wall_hi - λ[2] * d / 2) / bin_width[dim]) + 1

            K = CartesianIndex((NP[1:end-1]..., lo_idx))
            contr[1, dim] += rho_beads_K[K, mon_idx] * (exp(ε[1]) - 1.0)

            K = CartesianIndex((NP[1:end-1]..., hi_idx))
            contr[2, dim] += rho_beads_K[K, mon_idx] * (exp(ε[2]) - 1.0)       
        end
    end

    return contr
end