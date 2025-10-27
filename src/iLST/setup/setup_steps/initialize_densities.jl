
"""
    initialize_densities(molsys::MolecularSystem) -> BulkDensities

Constructs all relevant density fields for a given molecular system. These densities
serve as the initial state for bulk calculations in liquid-state theory or 
self-consistent field solvers.

This function:
- Converts species-level inputs (`pSpecies` or `density`) into normalized densities,
- Distributes total polymer density across all segments and accessible states,
- Computes segment–segment pair densities for all state combinations,
- Aggregates densities for individual monomers (beads),
- Tallies bond contributions based on segment connectivity and state pairs.

Arguments:
- `molsys`: A `MolecularSystem` containing sequences, state families, bonding
                   topology, and species-level thermodynamic inputs.

Returns:
- `BulkDensities`, with the following fields:
    - `species`   : Polymer density per species, normalized by segment count.
    - `segments`  : Matrix of state-resolved densities for each segment of each species.
    - `pairs`     : 4D array of segment–segment pair densities across all states.
    - `beads`     : Total density of each unique monomer (state-summed).
    - `bonds`     : Pairwise bond contributions separated into left/right beads.
"""



function initalize_densities(molsys :: MolecularSystem)

    @unpack configurations, properties = molsys
    @unpack monomers, species, bond_types = properties
    @unpack delta_muH = monomers
    @unpack species_charge = species

    n_species = length(species.species)

    
    rho_species  = zeros(n_species)
    rho_segments = Vector{Matrix{Float64}}(undef, n_species)
    rho_pairs    = Vector{Array{Float64, 4}}(undef, n_species)
    rho_beads    = zeros(length(species.monomers))
    rho_bonds    = zeros(2, length(bond_types))

    valences = get!(monomers, :valences, zeros(Float64, length(species.monomers)))

    @. species_charge = 0.0

    for (u, config) in enumerate(configurations)
        @unpack sequence, state_family = config

        n_segments = length(sequence)
        max_states = maximum(length(state_family[s]) for s in 1:n_segments)

        rho_segments[u] = zeros(max_states, n_segments)
        rho_pairs[u] = zeros(max_states, n_segments, max_states, n_segments)

        if !isnan(species.input_pSpecies[u])
            # Convert mol/L → reduced density
            conc = 10.0^(-species.input_pSpecies[u])
            density = concentration_to_density(conc, molsys.constants)
        else
            # Already given as reduced density
            density = species.input_densities[u]
        end

        polymer_density = density / n_segments

        rho_species[u] = polymer_density

        n_segments = length(sequence)

        for seg_i in 1:n_segments
            states_i = state_family[seg_i]
            Q_i = sum(@~ @. exp(-delta_muH[states_i]))

            # --- Update segment and bead densities ---
            for (idx_i, state_i) in enumerate(states_i)
                weight_i = exp(-delta_muH[state_i]) / Q_i
                rho_segments[u][idx_i, seg_i] += polymer_density * weight_i
                rho_beads[state_i] += polymer_density * weight_i
                species_charge[u] += weight_i * valences[state_i]

                
                # --- Pair and bond contributions ---
                for seg_j in 1:n_segments
                    states_j = state_family[seg_j]
                    Q_j = sum(@~ @. exp(-delta_muH[states_j]))

                    for (idx_j, state_j) in enumerate(states_j)
                        weight_j = exp(-delta_muH[state_j]) / Q_j
                        pair_density = polymer_density * weight_i * weight_j

                        # Update pair density
                        rho_pairs[u][idx_j, seg_j, idx_i, seg_i] += pair_density
                        
                        # Bond contribution if applicable
                        if isempty(config.topology.parents[seg_i])
                            # parent segment
                            continue
                        else
                            if config.topology.parents[seg_i][1] == seg_j
                                bond = (min(state_i, state_j), max(state_i, state_j))
                                idx = findfirst(==(bond), bond_types)
                                if idx !== nothing
                                    @. rho_bonds[:, idx] += pair_density
                                else
                                    error("Mismatch in bond from parent-child relation")
                                end
                            end
                        end
                    end
                end
            end
        end        
    end

    return BulkDensities(rho_species, rho_segments, rho_pairs, rho_beads, rho_bonds)
end