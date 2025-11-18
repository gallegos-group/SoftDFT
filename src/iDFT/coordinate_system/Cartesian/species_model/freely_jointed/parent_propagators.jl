function compute_parent_propagators(
        u,
        gC,
        bulk_system :: IsingLST, 
        geometry :: CartesianCoord, 
        fields :: SpatialFields, 
        species_model :: freely_jointed)


    # Related to species u
    config_u = bulk_system.molsys.configurations[u]

    # Fields
    @unpack lambda_K, lng_K, trapez = fields.excess

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

    gP = zeros(Float64, (NP..., max_states, max_states, n_segments )) # Parent propagator

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

    # Parent Propagator function
    for level in eachindex(levels)
        for (_, child) in enumerate(levels[level])
            for (idx_i, state_i) in enumerate(state_family[child])

                if fixed_segments[child]
                    K = CartesianIndex(fixed_coordinates[child])
                    gP[K, 1, idx_i, child] = 1.0
                else
                    for K in Rsys
                        gP[K, 1, idx_i, child] = 1.0
                    end
                end

                for parent in parents[child]
                    for (idx_j, state_j) in enumerate(state_family[parent])
                        bond, pair_i, pair_j = get_bond_info(state_i, state_j, bond_types)

                        Dij = (diameters[state_i] + diameters[state_j]) / 2.0
                        Dij2 = Dij^2
                        renorm = get_renorm(Dij, geometry)

                        @. f_hat = 0.0
                        if fixed_segments[parent]
                            K = CartesianIndex(fixed_coordinates[parent])
                            K_star = to_star_index(K, offset, mirrored)

                            calc_f_hat_parent!(f_hat, K, offset, mirrored, children, 
                                    parent, child, idx_j, state_j, pair_j, bond,
                                    lambda_K, lng_K, gP, gC)

                            coord = fixed_coordinates[parent]
                            for K1 in Rsys
                                dist2 = compute_distance_squared(K1, coord, bin_width)

                                if dist2 <= Dij2
                                    gP[K1, idx_j, idx_i, child] = abs(real(f_hat[K_star]))*exp(lng_K[K1, pair_i, bond]/2.0)/renorm/normative

                                    if isapprox(dist2, Dij2; atol = 1e-6)
                                        gP[K1, idx_j, idx_i, child] /= 2.0
                                    end
                                end
                            end
                        else
                            for K in Rsys
                                K_star = to_star_index(K, offset, mirrored)
                                calc_f_hat_parent!(f_hat, K, offset, mirrored, children, 
                                        parent, child, idx_j, state_j, pair_j, bond,
                                        lambda_K, lng_K, gP, gC)

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
                                gP[K, idx_j, idx_i, child] = abs(real(f_hat[K_star]))*exp(lng_K[K, pair_i, bond]/2.0)/normative
                            end
                        end
                    end
                end
            end

            for (segment, isfixed) in enumerate(fixed_segments)
                if isfixed
                    max_dist2 = species_model.max_distance[child, segment]^2
                    coord = fixed_coordinates[segment]
                    for K1 in Rsys
                        if compute_distance_squared(K1, coord, bin_width) > max_dist2
                            @. gP[K1, :, :, child] = 0.0
                        end
                    end
                end
            end
        end
    end

    return gP
end

function calc_f_hat_parent!(f_hat, K, offset, mirrored, children, parent, child, idx_j, state_j, pair_j, bond, lambda_K, lng_K, gP, gC)
    K_star = to_star_index(K, offset, mirrored)
    f_hat[K_star] = exp(-lambda_K[K, state_j] + lng_K[K, pair_j, bond]/2.0) *
               sum(@views gP[K, :, idx_j, parent])
    for (oc_index, other_child) in enumerate(children[parent])
        if other_child != child
            f_hat[K_star] *= sum(@views gC[K, :, idx_j, oc_index, parent])
        end
    end
end