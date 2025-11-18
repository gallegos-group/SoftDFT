function handle_fixed(
        u,
        gP,
        gC,
        bulk_system :: IsingLST, 
        geometry :: CartesianCoord, 
        fields :: SpatialFields, 
        species_model :: freely_jointed)

    # Fields
    @unpack lambda_K, trapez = fields.excess

    # Monomers
    @unpack delta_muH = bulk_system.molsys.properties.monomers

    # Geometry
    @unpack bin_width, NP = geometry

    # Topology
    @unpack state_family, topology = bulk_system.molsys.configurations[u]
    @unpack levels, children = topology

    # Fixed Segments
    fixed_density = fields.fixed[u].density
    fixed_segments = fields.fixed[u].segments
    fixed_coordinates = fields.fixed[u].coordinates

    n_segments = length(fixed_segments)

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
    
    if any(fixed_segments)
        norm = fixed_density[1]
        bulk_system.bulk.mu_species[u] = 0.0
        lambda = 0.0

        parent = findfirst(!, fixed_segments)
        @assert parent !== nothing "No non-fixed segment available for lambda normalization."

        for K in CartesianIndices(NP)
            for (idx_i, state_i) in enumerate(state_family[parent])

                contr = exp(-lambda_K[K, state_i])*sum(@views gP[K, :, idx_i, parent])

                if isempty(children[parent])
                    temp = sum(@views gC[K, :, idx_i, 1, parent])
                else
                    temp = 0.0
                    for (idx_c, child) in enumerate(children[parent])
                        temp1 = 0.0
                        for (idx_j, state_j) in enumerate(state_family[child])
                            temp2 = gC[K, idx_j, idx_c, idx_i, parent]
                            for (idx_oc, other_child) in enumerate(children[parent])
                                if other_child != child
                                    temp2 *= sum(@views gC[K, :, idx_oc, idx_i, parent])
                                end
                            end
                            temp1 += temp2
                        end
                        temp = temp1
                    end
                end
                lambda += contr*temp*trapez[K, state_i]*prod(bin_width)
            end
        end
        # println(lambda)
        # throw("")
        # parent = findfirst(fixed_segments)
        # lambda = 0.0
        # K = CartesianIndex(fixed_coordinates[1][parent])
        #     for (idx_i, state_i) in enumerate(state_family[parent])

        #         contr = exp(-mu_ex_K[K, state_i] - Ext[K, state_i] - delta_muH[state_i])*sum(@views gP[K, :, idx_i, parent])

        #         if isempty(children[parent])
        #             temp = sum(@views gC[K, :, idx_i, 1, parent])
        #         else
        #             temp = 0.0
        #             for (idx_c, child) in enumerate(children[parent])
        #                 temp1 = 0.0
        #                 for (idx_j, state_j) in enumerate(state_family[child])
        #                     temp2 = gC[K, idx_j, idx_c, idx_i, parent]
        #                     for (idx_oc, other_child) in enumerate(children[parent])
        #                         if other_child != child
        #                             temp2 *= sum(@views gC[K, :, idx_oc, idx_i, parent])
        #                         end
        #                     end
        #                     temp1 += temp2
        #                 end
        #                 temp = temp1
        #             end
        #         end
        #         lambda += contr*temp
        #     end
        #     println(lambda)
    else
        lambda = 1.0
        norm = (normative^(n_segments-1))
    end
    
    if !isfinite(lambda) || lambda <= 0
        error("Normalization error: λ = $lambda for species $u. Check propagator integrity.")
    end
    norm *= 1.0/lambda

    lambda /= normative^(n_segments-1)

    fields.fixed[u].lambda[1] = lambda

    return norm
end