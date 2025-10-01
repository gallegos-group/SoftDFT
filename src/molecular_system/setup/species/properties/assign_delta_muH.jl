"""
    assign_delta_muH!(
        properties::PropertiesStruct,
        used_monomers::Dict{String, Int},
        state_families::Dict{String, Vector{String}},
        constants::ConstantsStruct
    )

Computes and assigns the protonation free energy shifts Δμ_H for each monomer
in a state family based on reaction connectivity and chemical potentials.

### Description
For each pair of monomers in the same state family connected by a chemical reaction,
this function accumulates the free energy difference (Δμ) by:

- Starting from the pK of the reaction.
- Including contributions from context species (e.g., protons) via:
  - `pSpecies` if available
  - `-log(concentration)` if only density is known

The resulting Δμ values are assigned relative to the lexicographically first monomer in each state family,
which is chosen as the reference (Δμ = 0).

The output is stored in `properties.monomers[:delta_muH]`, indexed by `used_monomers`.

### Arguments
- `properties`: Contains monomer info (`:valences`, `:delta_muH`), species list, and input reactions.
- `used_monomers`: Maps each monomer name to its index in the monomer property arrays.
- `state_families`: Maps each monomer to its full equivalence class under protonation.
- `constants`: Includes unit conversions for computing concentrations from densities.

### Behavior
- Ignores any reaction not connecting two monomers from the same state family.
- Applies charge neutrality enforcement when valences are available.
- Automatically resolves context species contributions to Δμ from pK, pSpecies, or density.
- Modifies the reaction list in-place to remove equilibrium constraints already resolved.

### Errors
- Raises an error if a context species is not defined in the species list.
- Raises an error if a reaction involving state family members violates charge neutrality.

### Example
```julia
used = Dict("A" => 1, "B" => 2)
state_families = Dict("A" => ["A", "B"], "B" => ["A", "B"])
assign_delta_muH!(properties, used, state_families, constants)
"""


function assign_delta_muH!(
    properties::PropertiesStruct,
    used_monomers::Dict{String, Int},
    state_families::Dict{String, Vector{String}},
    constants::ConstantsStruct
)
    monomer_properties = properties.monomers
    species_data       = properties.species
    reactions          = properties.reactions

    log10 = log(10.0)
    delta_muH = zeros(length(used_monomers))

    idx_of    = used_monomers
    densities = species_data[:input_densities]
    pSpecies  = species_data[:input_pSpecies]
    names     = species_data[:species]
    name_to_index = Dict(name => i for (i, name) in enumerate(names))

    

    # Step 1: Determines if pSpecies can be defined
    #         for a neutralization reaction
    for rxn in reactions
        species, coeffs, pK = rxn.species, rxn.coeffs, rxn.pK[1]
        reactants = [species[i] for i in eachindex(coeffs) if coeffs[i] == -1]
        products  = [species[i] for i in eachindex(coeffs) if coeffs[i] == +1]

        # --- Neutralization case: ∅ → Y + Z or Y + Z → ∅ ---
        if isempty(reactants) || isempty(products)
            mon1, mon2 = isempty(reactants) ? products : reactants
            Δμ_total = -log10 * pK

            idx1, idx2 = name_to_index[mon1], name_to_index[mon2]
            ps1, ps2   = pSpecies[idx1], pSpecies[idx2]

            ps1_defined = !isnan(ps1)
            ps2_defined = !isnan(ps2)

            if ps1_defined
                δμ1 = -log10 * ps1
                δμ2 = Δμ_total - δμ1
                pSpecies[idx2] = δμ2 / -log10
            elseif ps2_defined
                δμ2 = -log10 * ps2
                δμ1 = Δμ_total - δμ2
                pSpecies[idx1] = δμ1 / -log10
            end
        end
    end


    # Step 2: Build Δμ graph from reactions (skips neutralization)
    cleaned_reactions = ReactionStruct[]
    Δμ_graph = Dict{String, Vector{Tuple{String, Float64}}}()
    for rxn in reactions
        species, coeffs, pK = rxn.species, rxn.coeffs, rxn.pK[1]
        reactants = [species[i] for i in eachindex(coeffs) if coeffs[i] == -1]
        products  = [species[i] for i in eachindex(coeffs) if coeffs[i] == +1]

        # === Charge neutrality check ===
        if haskey(monomer_properties, :valences)
            valences = monomer_properties[:valences]
            charge_change = 0.0
            for (i, sp) in enumerate(species)
                if haskey(used_monomers, sp)
                    idx = used_monomers[sp]
                    charge_change += coeffs[i] * valences[idx]
                end
            end
            if abs(charge_change) > 1e-6
                error("❌ Reaction $(species) violates charge neutrality: Δq = $charge_change")
            end
        end

        # Identify polymer monomers involved in a state family
        state_reactants = filter(s -> haskey(state_families, s), reactants)
        state_products = filter(s -> haskey(state_families, s), products)
        found_pair = false
        mon1, mon2 = "", ""
        fam1 = Vector{String}()

        for reactant in state_reactants
            for product in state_products
                mon1, mon2 = reactant, product
                fam1, fam2 = state_families[mon1], state_families[mon2]
                if fam1 == fam2
                    found_pair = true
                    break
                end
            end
            if found_pair
                break
            end
        end
        
        if !found_pair
            push!(cleaned_reactions, rxn)
            continue
        end
        
        # Determine direction: base → alt
        base = sort(fam1)[1]
        is_forward = (mon1 == base)
        from_mon, to_mon = is_forward ? (mon1, mon2) : (mon2, mon1)

        # Start from: Δμ = -log10 * pK
        Δμ = log10 * pK

        # Remove context species and accumulate their Δμ contribution
        new_species = String[]
        new_coeffs = Int[]
        for (i, s) in enumerate(species)
            if s == mon1 || s == mon2
                push!(new_species, s)
                push!(new_coeffs, coeffs[i])
            else
                idx = get(name_to_index, s, nothing)
                if isnothing(idx)
                    error("Context species $s is not defined in species list.")
                end
                coeff = coeffs[i]
                ps = pSpecies[idx]
                if !isnan(ps)
                    Δμ -= log10 * ps * coeff
                else
                    # We will need to self-constitently determine the correct Δμ
                    ρ = densities[idx]
                    if ρ > 0
                        c = density_to_concentration(ρ, constants)
                        Δμ += coeff * -log(c)
                    else
                        # Assuming pSpecies ≈ pK
                        Δμ = 0.0
                    end
                end
            end
        end
        
        push!(get!(Δμ_graph, from_mon, Vector{Tuple{String, Float64}}()), (to_mon, Δμ))
        push!(get!(Δμ_graph, to_mon, Vector{Tuple{String, Float64}}()), (from_mon, Δμ))

        # DO not need to keep this reaction because chemical equilibrium is automatically satisfied
        @. rxn.coeffs = 0
        # push!(cleaned_reactions, ReactionStruct(species, fill!(coeffs, 0.0), [pK]))
    end

    # Step 2: Traverse each state family and accumulate Δμ from base
    for (_, family) in state_families
        base = sort(family)[1]
        visited = Set{String}()
        work = Dict(base => 0.0)

        function dfs(monomer, acc)
            push!(visited, monomer)
            for (neighbor, dμ) in get(Δμ_graph, monomer, [])
                if neighbor ∉ visited
                    work[neighbor] = acc + dμ
                    dfs(neighbor, acc + dμ)
                end
            end
        end
        

        dfs(base, 0.0)
        
        for (mon, val) in work
            delta_muH[idx_of[mon]] = val
        end
        delta_muH[idx_of[base]] = 0.0
    end
    
    # Step 1: Determines if pSpecies can be defined
    #         for a neutralization reaction
    for rxn in reactions
        species, coeffs, pK = rxn.species, rxn.coeffs, rxn.pK[1]
        reactants = [species[i] for i in eachindex(coeffs) if coeffs[i] == -1]
        products  = [species[i] for i in eachindex(coeffs) if coeffs[i] == +1]

        # --- Neutralization case: ∅ → Y + Z or Y + Z → ∅ ---
        if (isempty(reactants) || isempty(products)) && all(c != 0 for c in coeffs)
            mon1, mon2 = isempty(reactants) ? products : reactants

            idx1, idx2 = name_to_index[mon1], name_to_index[mon2]
            ps1, ps2   = pSpecies[idx1], pSpecies[idx2]

            ps1_defined = !isnan(ps1)
            ps2_defined = !isnan(ps2)

            if ps1_defined || ps2_defined
                @. rxn.coeffs = 0
            end
        end
    end

    monomer_properties[:delta_muH] = delta_muH
    # empty!(reactions)
    # append!(reactions, cleaned_reactions)
end