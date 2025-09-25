"""
    parse_reactions(dataset::Dict) -> Vector{ReactionStruct}

Parses chemical reactions from the YAML input dataset.

Each reaction entry must contain:
- `"species"`: a list of species names involved in the reaction (e.g., `["A", "H", "B"]`)
- `"coeffs"`: stoichiometric coefficients corresponding to each species (e.g., `[-1, -1, +1]`)
- `"pK"`: a single equilibrium constant for the reaction, given in base-10 logarithmic form

**Constraints:**
- The number of species must match the number of coefficients.
- All coefficients must be exactly ±1, either as integers or floating-point values (e.g., `1`, `-1`, `1.0`, `-1.0`).
- Coefficients like `0`, `2`, `0.5`, or `-2` are not allowed.

Returns:
- A vector of `ReactionStruct` structs with validated and normalized coefficients.

Raises:
- An error if any reaction is ill-formed, uses invalid coefficients, or violates supported reaction formats.

Example:
```yaml
reactions:
  - species: ["A", "H", "B"]
    coeffs: [-1, -1, 1]
    pK: 5.0

"""


function parse_reactions(dataset::Dict)
    if !haskey(dataset, "reactions") || isnothing(dataset["reactions"])
        return ReactionStruct[]
    end

    raw_reactions = dataset["reactions"]

    if !(raw_reactions isa Vector)
        error("Expected 'reactions' to be a list of reaction entries, but got: $(typeof(raw_reactions))")
    end
    
    reaction_entries = ReactionStruct[]

    for rxn in raw_reactions
        species = rxn["species"]
        coeffs_raw = rxn["coeffs"]
        pK = rxn["pK"]

        if length(species) != length(coeffs_raw)
            error("❌ Mismatch between number of species and coefficients in reaction: $rxn")
        end

        coeffs = Int[]
        for c in coeffs_raw
            # Convert to Int if it's a Float like 1.0 or -1.0
            c_int = round(Int, c)
            if abs(c - c_int) > 1e-8
                error("❌ Non-integer stoichiometric coefficient $c in reaction: $rxn. Only ±1 are supported.")
            end
            if !(c_int == 1 || c_int == -1)
                error("❌ Unsupported coefficient $c_int in reaction: $rxn. Only ±1 are allowed.")
            end
            push!(coeffs, c_int)
        end

        push!(reaction_entries, ReactionStruct(species, coeffs, [pK]))
    end

    return reaction_entries
end


