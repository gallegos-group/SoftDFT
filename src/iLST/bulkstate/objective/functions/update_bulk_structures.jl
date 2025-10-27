"""
    update_bulk_structures!(molsys, bulk, solution)

Updates `rho.species` in `bulk` using the current solver values of `log(ρ_i)`.

- The input `solution` is a tuple: (log_rho_rxn, log_rho_pH)
- For each adjustable species, reconstructs ρ_i = exp(log(ρ_i))
- Fixed-density species are not updated.
"""
function update_bulk_structures!(molsys::MolecularSystem, bulk::BulkState, solution)
    @unpack species, input_pSpecies, adjustable = molsys.properties.species
    @unpack rho = bulk
    reactions = molsys.properties.reactions

    log_rho_rxn, log_rho_pH, Q_solver = solution

    # Set Q_solver to control charge-balancing species
    balance_charge!(molsys, bulk, Q_solver[1])

    # Identify species involved in meaningful reactions
    species_in_reactions = Set{String}()
    for rxn in reactions
        if all(c == 0 for c in rxn.coeffs)
            continue  # bookkeeping reaction
        end
        for s in rxn.species
            push!(species_in_reactions, s)
        end
    end

    idx_rxn = 1
    idx_pH  = 1

    for (i, s) in enumerate(species)
        if !isnan(input_pSpecies[i])
            rho.species[i] = exp(log_rho_pH[idx_pH])
            idx_pH += 1
        elseif s in species_in_reactions
            rho.species[i] = exp(log_rho_rxn[idx_rxn])
            idx_rxn += 1
        end
    end
end

