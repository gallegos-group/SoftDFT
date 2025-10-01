"""
    build_used_monomers(species_data, reaction_entries) -> Dict{String, Int}

Constructs a dictionary mapping each relevant monomer name to a unique integer index.

### Monomer Inclusion Criteria:
- All single-letter monomers that appear in the `sequence` field of each species (typically from polymers).
- All species referenced in any chemical reaction, regardless of whether they are used in a sequence.
  This allows for full construction of state families in a later step.

### Behavior:
- Reactions are NOT filtered.
- Returns a dictionary of monomer names sorted alphabetically with 1-based indices.
- The dictionary is used to index monomer-specific properties like valence, diameter, and Δμ_H.

### Arguments
- `species_data::Vector{Dict}`: Parsed species section from the YAML input.
- `reaction_entries::Vector{ReactionStruct}`: Parsed list of chemical reactions. These are not mutated.

### Returns
- `Dict{String, Int}`: Mapping from monomer name to unique index.

### Example
If `"sequence": "AAA"` and a reaction `["A", "Y", "B"]` exist, the result might be:
    Dict("A" => 1, "B" => 2, "Y" => 3)
"""
function build_used_monomers(species_data::Vector{Dict{Any, Any}}, reaction_entries::Vector{ReactionStruct})
    monomer_names = Set{String}()

    # Step 1: Extract monomers from species sequences
    for entry in species_data
        seq = entry["sequence"]
        for c in seq
            if isletter(c)
                push!(monomer_names, string(c))
            end
        end
    end

    # Step 2: Add all species that appear in any reaction
    for rxn in reaction_entries
        for s in rxn.species
            push!(monomer_names, s)
        end
    end

    # Step 3: Assign sorted indices
    sorted_names = sort(collect(monomer_names))
    return Dict(name => i for (i, name) in enumerate(sorted_names))
end