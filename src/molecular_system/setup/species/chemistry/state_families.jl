function build_state_families(species_entries::Vector{SpeciesEntry}, reactions::Vector{ReactionStruct})

    # Step 1b: Assign neutralizer roles based on reactions
    assign_neutralizer_roles!(species_entries, reactions)

    # Step 1: Extract roles
    species_roles = Dict{String, Symbol}()
    for entry in species_entries
        species_roles[entry.name] = entry.role
    end

    # Step 2: Initialize connectivity and classification sets
    neighbors = Dict{String, Set{String}}()
    pending_edges = Vector{Tuple{Vector{String}, Vector{String}, ReactionStruct}}()
    neutralizers = Set{String}(filter(x -> species_roles[x] == :neutralizer, keys(species_roles)))
    partners      = Set{String}(filter(x -> species_roles[x] == :partner, keys(species_roles)))
    monomers      = Set{String}(filter(x -> species_roles[x] == :monomer, keys(species_roles)))
    polymer_monomers = find_polymer_monomers(species_entries)

    # Step 3: Parse reactions
    for rxn in reactions
        species, coeffs = rxn.species, rxn.coeffs
        reactants = [species[i] for i in eachindex(coeffs) if coeffs[i] == -1]
        products  = [species[i] for i in eachindex(coeffs) if coeffs[i] == +1]

        # Skip neutralization reactions (∅ → X or X → ∅)
        if isempty(reactants) || isempty(products)
            for s in (isempty(reactants) ? products : reactants)
                if s in polymer_monomers
                    error("❌ Monomer '$s' appears in a neutralization reaction, but it is part of a polymer sequence. This is not allowed.")
                end
            end
            continue
        end


        # Validate supported reaction form
        if !is_supported_reaction_form(reactants, products)
            error("❌ Unsupported reaction form: $(rxn.species) with coeffs $(rxn.coeffs)")
        end

        # For reactions of the form A + Y -> B
        if length(reactants) == 2 && length(products) == 1
            roles = get(species_roles, reactants[1], :monomer), get(species_roles, reactants[2], :monomer)
            if all(r -> r == :monomer, roles)
                error("❌ Reaction $(rxn.species) has two reactants and one product, but both reactants are marked or assumed as :monomer. Please declare one (e.g. '$(reactants[2])') as a :partner or :neutralizer.")
            end
        end

        # For reactions of the form A -> B + Y
        if length(reactants) == 1 && length(products) == 2
            roles = get(species_roles, products[1], :monomer), get(species_roles, products[2], :monomer)
            if all(r -> r == :monomer, roles)
                error("❌ Reaction $(rxn.species) has one reactant and two products, but both products are marked or assumed as :monomer. Please declare one (e.g. '$(products[2])') as a :partner or :neutralizer.")
            end
        end


        # Extract only monomers for connectivity
        lhs_candidates = filter(s -> get(species_roles, s, :monomer) ∉ (:partner, :neutralizer), reactants)
        rhs_candidates = filter(s -> get(species_roles, s, :monomer) ∉ (:partner, :neutralizer), products)

        total_candidates = length(lhs_candidates) + length(rhs_candidates)

        if total_candidates > 2
            error("❌ Reaction $(rxn.species) connects more than two candidate monomers — ambiguous state family inference.")
        elseif total_candidates == 2
            a = lhs_candidates[1]
            b = rhs_candidates[1]
            neighbors[a] = get(neighbors, a, Set{String}())
            neighbors[b] = get(neighbors, b, Set{String}())
            push!(neighbors[a], b)
            push!(neighbors[b], a)
        elseif total_candidates == 1
            # Ambiguous — could be monomer + partner → unknown product
            push!(pending_edges, (reactants, products, rxn))
        else
            # No valid monomer candidates — ignore
            continue
        end
    end

    # Step 4: Resolve ambiguous cases
    resolve_pending_edges!(pending_edges, neighbors, species_roles)

    # Step 5: Finalize state families
    all_monomers = find_all_monomers(species_entries)

    return finalize_state_families(neighbors, all_monomers)
end


function resolve_pending_edges!(
    pending_edges::Vector{Tuple{Vector{String}, Vector{String}, ReactionStruct}},
    neighbors::Dict{String, Set{String}},
    species_roles::Dict{String, Symbol}
)
    for (reactants, products, rxn) in pending_edges
        candidates = filter(s -> species_roles[s] == :monomer, reactants ∪ products)

        if length(candidates) != 2
            continue  # skip: not enough monomers to connect
        end

        a, b = candidates
        fam_a = find_component_membership(neighbors, a)
        fam_b = find_component_membership(neighbors, b)

        if !haskey(neighbors, a)
            neighbors[a] = Set{String}()
        end
        if !haskey(neighbors, b)
            neighbors[b] = Set{String}()
        end

        if a in fam_b || b in fam_a
            continue  # already connected
        else
            push!(neighbors[a], b)
            push!(neighbors[b], a)
        end
    end
end

function find_all_monomers(species_entries::Vector{SpeciesEntry})
    monomers = Set{String}()
    for entry in species_entries
        seq = entry.name  # assuming .name is the sequence string like "CCC"
        for c in seq
            if isletter(c)
                push!(monomers, string(c))
            end
        end
    end
    return monomers
end

"""
    finalize_state_families(
        neighbors::Dict{String, Set{String}},
        all_monomers::Set{String}
    ) -> Dict{String, Vector{String}}

Builds the final mapping of monomers to their sorted state families.

Each state family is a connected component in the undirected graph `neighbors`.
The returned dictionary maps each monomer to its full state family (as a sorted `Vector{String}`).

### Arguments:
- `neighbors`: Graph of pairwise connections built from resolved reactions
- `all_monomers`: Set of all monomers declared in the `species:` block of the YAML input

### Returns:
- `Dict{String, Vector{String}}` where each monomer is mapped to its state family

Monomers not connected to any others will form singleton state families (i.e., `"Y" => ["Y"]`).
"""


function finalize_state_families(neighbors::Dict{String, Set{String}}, all_monomers::Set{String})
    visited = Set{String}()
    state_families = Dict{String, Vector{String}}()

    function dfs(root::String)
        stack = [root]
        group = String[]

        while !isempty(stack)
            node = pop!(stack)
            if node ∉ visited
                push!(visited, node)
                push!(group, node)
                append!(stack, collect(get(neighbors, node, Set())))
            end
        end

        return sort(group)
    end

    for mon in all_monomers
        if mon ∉ visited
            group = dfs(mon)
            for m in group
                state_families[m] = group
            end
        end
    end

    return state_families
end


"""
    is_supported_reaction_form(reactants::Vector{String}, products::Vector{String}) -> Bool

Returns `true` if the given reaction is of a supported form for state family inference.

### Supported reaction forms:
- `A → B` (1 → 1)
- `A → B + Y` (1 → 2)
- `A + Y → B` (2 → 1)

Where `Y` may be a context species (e.g., proton, water, etc.).

### Returns:
- `true` if the reaction is in one of the supported (1,1), (1,2), or (2,1) forms.
- `false` otherwise (e.g., 2 → 2, 3 → 1, or 1 → 3 are invalid).

This check is used to restrict which reactions can be interpreted as reversible state transitions for building monomer state families.
"""

function is_supported_reaction_form(reactants::Vector{String}, products::Vector{String})
    r, p = length(reactants), length(products)
    return (r == 1 && p in (1, 2)) || (r == 2 && p == 1)
end


function assign_neutralizer_roles!(
    species_entries::Vector{SpeciesEntry},
    reactions::Vector{ReactionStruct}
)
    # Detect species that participate in ∅ → X + Y or X + Y → ∅
    neutralizer_names = detect_neutralization_participants(reactions)

    for entry in species_entries
        if entry.name in neutralizer_names
            entry.role = :neutralizer
        end
    end
end



"""
    detect_neutralization_participants(reactions) -> Set{String}

Identifies monomers involved in neutralization reactions, defined as
having either an empty reactant or product side.

Returns a set of neutralizer species.
"""
function detect_neutralization_participants(reactions::Vector{ReactionStruct})
    neutralizers = Set{String}()

    for rxn in reactions
        species, coeffs = rxn.species, rxn.coeffs
        reactants = [species[i] for i in eachindex(coeffs) if coeffs[i] == -1]
        products  = [species[i] for i in eachindex(coeffs) if coeffs[i] == +1]

        if isempty(reactants) || isempty(products)
            for s in (isempty(reactants) ? products : reactants)
                push!(neutralizers, s)
            end
        end
    end

    return neutralizers
end


"""
    find_polymer_monomers(species_entries) -> Set{String}

Returns a set of monomer names that appear in polymer sequences of length > 1.
These monomers are considered part of polymers.
"""
function find_polymer_monomers(species_entries::Vector{SpeciesEntry})
    monomers = Set{String}()
    for entry in species_entries
        seq = entry.name  # species name is the sequence
        if length(seq) > 1
            for c in seq
                if isletter(c)
                    push!(monomers, string(c))
                end
            end
        end
    end
    
    return monomers
end

