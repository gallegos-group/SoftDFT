function compute_child_propagators(
        u, 
        bulk_system :: IsingLST, 
        geometry :: CartesianCoord, 
        fields :: SpatialFields, 
        species_model :: freely_jointed)

    # Related to species u
    config_u = bulk_system.molsys.configurations[u]

    # Fields
    @unpack mu_ex_K, lng_K, Ext, trapez = fields.excess

    # Monomers
    @unpack diameters, delta_muH = bulk_system.molsys.properties.monomers

    # FFT
    @unpack f_hat, plan_forward, plan_backward, weight_bonds_hat = fields.fourier

    # Geometry
    @unpack NP, bin_width, mirrored, offset = geometry
    Rsys = CartesianIndices(NP)

    # Topology
    @unpack sequence, state_family, topology = config_u
    @unpack parents, children, levels = topology
    @unpack bond_types = bulk_system.molsys.properties

    n_segments = length(sequence)

    # Fixed Segments
    fixed_species_u = fields.fixed[u]
    fixed_segments = fixed_species_u.segments
    fixed_coordinates = fixed_species_u.coordinates
    
    max_states = maximum(length, state_family)
    max_childs = max(maximum(x -> length(x), children), 1)

    gC = zeros(Float64, (NP..., max_states, max_childs, max_states, n_segments)) # Child propagator

    normative = 1.0
    for level in reverse(eachindex(levels))
        for (_, parent) in enumerate(levels[level])
            temp = 0.0
            for (idx_i, state_i) in enumerate(state_family[parent])
                temp += exp(-delta_muH[state_i])
            end
            normative *= temp^(1.0/n_segments)
        end
    end

    # Child propagator functions
    for level in reverse(eachindex(levels))
        for (_, parent) in enumerate(levels[level])
            for (idx_i, state_i) in enumerate(state_family[parent])
                if fixed_segments[parent]
                    K = CartesianIndex(fixed_coordinates[parent])
                    gC[K, 1, 1, idx_i, parent] = 1.0
                else
                    for K in Rsys
                        gC[K, 1, 1, idx_i, parent] = 1.0
                    end
                end
        
                for (idx_c, child) in enumerate(children[parent])
                    for (idx_j, state_j) in enumerate(state_family[child])
                        bond, pair_i, pair_j = get_bond_info(state_i, state_j, bond_types)

                        Dij = (diameters[state_i] + diameters[state_j]) / 2.0
                        renorm = get_renorm(Dij, geometry)
                        Dij2 = Dij^2                      

                        @. f_hat = 0.0
                        if fixed_segments[parent]
                            coord = fixed_coordinates[parent]
                            K = CartesianIndex(coord)

                            if fixed_segments[child]
                                K_child = CartesianIndex(fixed_coordinates[child])
                                K_star_child = to_star_index(K_child, offset, mirrored)
                                calc_f_hat_child!(f_hat, K_child, offset, mirrored, children, child, idx_j, state_j, pair_j, bond, mu_ex_K, Ext, lng_K, gC)
                                gC[K, idx_j, idx_c, idx_i, parent] = abs(real(f_hat[K_star_child]))*exp(lng_K[K, pair_i, bond]/2.0 - delta_muH[state_j])/renorm/prod(bin_width)/normative
                            else    
                                gC[K, idx_j, idx_c, idx_i, parent] = 0.0
                                for K1 in Rsys
                                    dist2 = compute_distance_squared(K1, coord, bin_width)

                                    if dist2 <= Dij2
                                        K1_star = to_star_index(K1, offset, mirrored)
                                        calc_f_hat_child!(f_hat, K1, offset, mirrored, children, child, idx_j, state_j, pair_j, bond, mu_ex_K, Ext, lng_K, gC)

                                        temp = abs(real(f_hat[K1_star]))*exp(lng_K[K, pair_i, bond]/2.0 - delta_muH[state_j])*prod(bin_width)*trapez[K1, state_j]/renorm

                                        if isapprox(dist2, Dij2; atol=1e-6)
                                            temp /= 2.0
                                        end

                                        gC[K, idx_j, idx_c, idx_i, parent] += temp/normative
                                    end
                                end
                            end
                        elseif fixed_segments[child]
                            K = CartesianIndex(fixed_coordinates[child])
                            K_star = to_star_index(K, offset, mirrored)
                            calc_f_hat_child!(f_hat, K, offset, mirrored, children, child, idx_j, state_j, pair_j, bond, mu_ex_K, Ext, lng_K, gC)
                            coord = fixed_coordinates[child]

                            for K1 in Rsys
                                dist2 = compute_distance_squared(K1, coord, bin_width)

                                if dist2 <= Dij2
                                    gC[K1, idx_j, idx_c, idx_i, parent] = abs(real(f_hat[K_star]))*exp(lng_K[K1, pair_i, bond]/2.0 - delta_muH[state_j])/renorm/normative
                                    if isapprox(dist2, Dij2; atol=1e-6)
                                        gC[K1, idx_j, idx_c, idx_i, parent] /= 2.0
                                    end
                                end
                            end
                        else
                            for K in Rsys
                                K_star = to_star_index(K, offset, mirrored)
                                calc_f_hat_child!(f_hat, K, offset, mirrored, children, child, idx_j, state_j, pair_j, bond, mu_ex_K, Ext, lng_K, gC)
                                
                                f_hat[K_star] *= trapez[K, state_j]
                            end

                            apply_mirroring!(f_hat, geometry)   

                            f_hat = plan_forward * f_hat

                            apply_filter!(f_hat, geometry)
                            
                            for K_star in CartesianIndices(f_hat)
                                f_hat[K_star] *= weight_bonds_hat[K_star, bond]
                            end

                            f_hat = plan_backward * f_hat

                            for K in Rsys
                                K_star = to_star_index(K, offset, mirrored)
                                gC[K, idx_j, idx_c, idx_i, parent] = abs(real(f_hat[K_star]))*exp(lng_K[K, pair_i, bond]/2.0 - delta_muH[state_j])/normative
                            end
                        end
                    end
                end
            end

            for (segment, isfixed) in enumerate(fixed_segments)
                if isfixed
                    max_dist2 = species_model.max_distance[parent, segment]^2
                    coord = fixed_coordinates[segment]
                    for K1 in Rsys
                        if compute_distance_squared(K1, coord, bin_width) > max_dist2
                           @. gC[K1, :, :, :, parent] = 0.0
                        end
                    end
                end
            end
        end
    end

    return gC
end


@inline function compute_distance_squared(K::CartesianIndex, coord::NTuple{N, Int}, bin_width) where {N}
    idx = Tuple(K)
    d2 = 0.0
    @inbounds @simd for v in 1:N
        Δ = (idx[v] - coord[v]) * bin_width[v]
        d2 += Δ^2
    end
    return d2
end

@inline function compute_distance_squared(K::CartesianIndex, coord::CartesianIndex, bin_width)
    idx = Tuple(K)
    idx1 = Tuple(coord)
    d2 = 0.0
    @inbounds @simd for v in eachindex(idx)
        Δ = (idx[v] - idx1[v]) * bin_width[v]
        d2 += Δ^2
    end
    return d2
end


function calc_f_hat_child!(f_hat, K, offset, mirrored, children, child, idx_j, state_j, pair_j, bond, mu_ex_K, Ext, lng_K, gC)
    K_star = to_star_index(K, offset, mirrored)
    f_hat[K_star] = exp(-mu_ex_K[K, state_j] - Ext[K, state_j] + lng_K[K, pair_j, bond] / 2.0)
    for (idx_gc, grandchild) in enumerate(children[child])
        f_hat[K_star] *= sum(@views gC[K, :, idx_gc, idx_j, child])
    end
end

"""
    get_bond_info(state_i, state_j, bond_types) -> bond_index, pair_i, pair_j

Returns:
- `bond_index`: index of the (ordered) state pair in `bond_types`
- `pair_i`, `pair_j`: indices (1 or 2) indicating where `state_i` and `state_j` appear in the ordered pair

If the pair is not found in `bond_types`, returns `nothing` for `bond_index`.
"""
function get_bond_info(state_i, state_j, bond_types)
    if state_i <= state_j
        pair = (state_i, state_j)
        pair_i, pair_j = 1, 2
    else
        pair = (state_j, state_i)
        pair_i, pair_j = 2, 1
    end

    bond_index = findfirst(==(pair), bond_types)
    return bond_index, pair_i, pair_j
end


function get_renorm(Dij, geometry)
    if geometry isa CartesianZ
        return 2.0*Dij
    elseif geometry isa CartesianXYZ
        return 4.0*pi*Dij^2
    else
        error("Not added yet.")
    end
end