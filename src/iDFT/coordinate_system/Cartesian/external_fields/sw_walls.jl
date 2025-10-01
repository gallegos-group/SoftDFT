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

function eval_External(molsys, geometry::CartesianCoord, fields, ::Val{:sw_walls})
    @unpack trapez, Ext = fields.excess
    @unpack NP, bin_width, features = geometry
    @unpack mirrored, offset = features
    @unpack diameters = molsys.properties.monomers

    sw_walls = features[:external_field][:sw_walls]
    energys = sw_walls[:surface_energys]
    lambdas = sw_walls[:surface_lambdas]
    wall_positions = sw_walls["position"]

    Rsys = CartesianIndices(NP)

    for K in Rsys
        idx = Tuple(K)

        for (dim, walls) in enumerate(wall_positions)
            isempty(walls) && continue

            wall_lo = walls[1]
            wall_hi = walls[2]

            for mon_idx in eachindex(diameters)
                d = diameters[mon_idx]
                λ = lambdas[dim][mon_idx]
                ε = energys[dim][mon_idx]

                
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


function contact_value_theorem(molsys, geometry::CartesianCoord, fields, ::Val{:sw_walls})
    @unpack Ext = fields.excess
    @unpack NP, bin_width, features = geometry
    @unpack mirrored, offset = features
    @unpack diameters = molsys.properties.monomers

    sw_walls = features[:external_field][:sw_walls]
    energys = sw_walls[:surface_energys]
    lambdas = sw_walls[:surface_lambdas]
    wall_positions = sw_walls["position"]

    rho_beads_K = fields.rho_K.beads
    
    contr = zeros(Float64, (2,length(geometry.NP)))
    for (dim, walls) in enumerate(wall_positions)
        isempty(walls) && continue

        wall_lo = walls[1]
        wall_hi = walls[2]

        for mon_idx in eachindex(diameters)
            d = diameters[mon_idx]
            λ = lambdas[dim][mon_idx]
            ε = energys[dim][mon_idx]

            
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