function compute_neutrality_shift(
        bulk_system :: IsingLST, 
        geometry :: CartesianCoord, 
        fields :: SpatialFields)

    @unpack configurations = bulk_system.molsys
    @unpack valences = bulk_system.molsys.properties.monomers
    @unpack NP, bin_width, total_charge = geometry
    @unpack fixed = fields
    @unpack trapez = fields.excess

    rho_segments_K = fields.rho_K.segments
    bin = prod(bin_width)
    
    Rsys = CartesianIndices(NP)

    Qtot = sum(total_charge)

    Num_species = zeros(Float64, length(configurations))
    species_charge = zeros(Float64, length(configurations))
    for (u, config) in enumerate(configurations)
        rho_segments_u = rho_segments_K[u]
        for seg_i in eachindex(config.sequence)
            for (idx_i, state_i) in enumerate(config.state_family[seg_i])
                if fixed[u].segments[seg_i]
                    K = CartesianIndex(fixed[u].coordinates[1][seg_i])
                    num = rho_segments_u[K, idx_i, seg_i]
                else
                    num = 0.0
                    for K in Rsys
                        num += rho_segments_u[K, idx_i, seg_i]*trapez[K, state_i]*bin
                    end
                end

                Num_species[u] += num
                species_charge[u] += num*valences[state_i]
            end
        end

        Num_species[u] /= length(config.sequence)
        
        if Num_species[u] <= 0.0
            species_charge[u] = 0.0
        else
            species_charge[u] /= Num_species[u]
        end
    end

    for (u, config) in enumerate(configurations)
        has_state_family = any(length(fam) > 1 for fam in config.state_family)
        has_fixed_segment = any(fields.fixed[u].segments)
        isevaluated = bulk_system.bulk.evaluation[u] isa SimulationEval

        if has_state_family || has_fixed_segment || isevaluated
            Qtot += Num_species[u] * species_charge[u]
            species_charge[u] = 0.0
        end
    end
    
    err = abs(Qtot + sum(@~ @. Num_species * species_charge))
err1 = err
    tolerance = 1e-10
    max_iters = 1000
    Psi1C = 0.0
    iter = 0
    Cstar = exp(-Psi1C)

    iter = 1
    while err > tolerance
        fx = Qtot
        fx1 = 0.0
        for i in eachindex(species_charge)
            t = (Cstar ^ species_charge[i]) * Num_species[i] * species_charge[i]
            fx += t
            fx1 -= t * species_charge[i]
        end

        err = abs(fx)

        # Ensure fx1 is not too small to prevent division by zero
        if abs(fx1) < 1e-10
            # println("Derivative too small, stopping iteration.")
            break
        end

        del = max(1.0, abs(fx / fx1))
        Psi1C -= fx / (fx1 * del)

        Cstar = exp(-Psi1C)

        iter += 1
        if iter > max_iters
            # println("Failed to converge in potential after $max_iters iterations.")
            break
        end
    end
    
    return @. fields.excess.PsiC + Psi1C
end