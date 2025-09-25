

"""
    balance_charge!(molsys, bulk; strategy = :equal)

Enforces electroneutrality by adjusting densities of species marked as `compensator = true`.

The total charge introduced by the solver is redistributed across those compensators,
respecting lower bounds defined by `input_densities`, and updating `adjustable` flags
accordingly.

Arguments:
- `strategy`: `:equal` (default) or `:valence_weighted` distribution of charge

Modifies:
- `rho.species`
- `adjustable`
"""


function balance_charge!(
    molsys::MolecularSystem,
    bulk::BulkState,
    Q_solver;
    strategy::Symbol = :equal  # or :valence_weighted
)

    @unpack rho = bulk
    @unpack species, input_densities, compensator, 
            adjustable, species_charge = molsys.properties.species

    # Step 2: Decide which compensators absorb this charge
    if iszero(Q_solver)
        # Choose one "side" arbitrarily — pick negative compensators
        sign_Q = -1  # or +1, but must be consistent

        @. adjustable = false
        for i in eachindex(rho.species)
            if compensator[i] && sign(species_charge[i]) == -sign_Q
                adjustable[i] = true
            end
        end
        
        return
    end

    # Otherwise: normal logic when Q_solver ≠ 0
    sign_Q = sign(Q_solver)
    @. adjustable = false
    for i in eachindex(rho.species)
        if compensator[i] && sign(species_charge[i]) == -sign_Q
            adjustable[i] = true
        end
    end
    
    if !any(adjustable)
        error("⚠️ No compensating species available to absorb net charge $Q_solver")
    end
    

    # Step 3: Compute total weight over adjustable species
    total_weight = 0.0
    for i in eachindex(rho.species)
        if adjustable[i]
            total_weight += strategy == :valence_weighted ? abs(species_charge[i]) : 1.0
        end
    end

    # Step 4: Redistribute charge directly and update densities
    for i in eachindex(rho.species)
        if adjustable[i]
            weight = strategy == :valence_weighted ? abs(species_charge[i]) : 1.0
            Δ = -Q_solver * (weight / total_weight) / species_charge[i]
            rho.species[i] = input_densities[i] + max(Δ, 0.0)
        elseif compensator[i]
            rho.species[i] = input_densities[i]
        end
    end
end
