"""
    parse_species(dataset::Dict) -> Vector{SpeciesEntry}

Parses the `species` section of the input dataset. Each entry must define:

- `"sequence"`: the name or sequence identifier of the species (required)

Optional fields:
- `"density"`: fixed number density in reduced (simulation) units
- `"concentration"`: fixed concentration in mol/L (Molar units)
- `"pSpecies"`: Sets the chemical potential via μ = ln(10) × pSpecies
- `"adjustable"`: marks the species as a degree of freedom to satisfy electroneutrality

**Rules enforced:**
- A species may define at most one of: `density`, `concentration`, or `pSpecies`.
- It is an error to define both `density` and `concentration`.
- It is an error to define both `pSpecies` and `adjustable`.

**Note:**
Species with no specified `density`, `concentration`, `pSpecies`, or `adjustable` are allowed.
These under-specified species may be involved in reactions and resolved later.

Returns a vector of `SpeciesEntry` structs, one per species.

Throws an error if:
- The `species` section is missing
- An entry lacks the required `sequence` field
- An entry violates any exclusivity constraint
"""

mutable struct SpeciesEntry
    name::String
    density::Float64
    concentration :: Float64
    pSpecies::Float64   # Use NaN if not present
    adjustable::Bool
    role::Symbol
end

function parse_species(dataset::Dict)::Vector{SpeciesEntry}
    if !haskey(dataset, "species")
        error("Input file is missing the required 'species' section.")
    end

    species_entries = SpeciesEntry[]

    for entry in dataset["species"]
        has_sequence = haskey(entry, "sequence")
        has_density  = haskey(entry, "density")
        has_concentration  = haskey(entry, "concentration")
        has_pSpecies = haskey(entry, "pSpecies")
        adjustable = get(entry, "adjustable", false)
        role = Symbol(get(entry, "role", "monomer"))

        if has_sequence
            name = entry["sequence"]
        else
            error("❌ Each species must define a 'sequence' field.")
        end

        # Disallow both being defined
        if (has_density || has_concentration) && has_pSpecies
            error("❌ Species '$name' is overspecified: both 'density/concentration' and 'pSpecies' are defined. Only one may be specified.")
        elseif has_density && has_concentration
            error("❌ Species '$name' is overspecified: both 'density' and 'concentration' Only one may be specified.")
        end
    
        # Disallow pSpecies + adjustable
        if has_pSpecies && adjustable
            error("❌ Species '$name' has 'pSpecies' and is also marked as 'adjustable'. This is not allowed.")
        end
        
        # Use fallback values
        density   = get(entry, "density", 0.0)
        concentration   = get(entry, "concentration", 0.0)
        pSpecies  = get(entry, "pSpecies", NaN)
    
        push!(species_entries, SpeciesEntry(name, density, concentration, pSpecies, adjustable, role))
    end
    

    return species_entries
end