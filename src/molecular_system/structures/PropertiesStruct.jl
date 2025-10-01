
"""
    PropertiesStruct

Stores all physical, thermodynamic, and interaction parameters for the system.

Fields:
- `monomers`    : Dict of monomer-level constants (e.g., "diameters", "valences", "delta_muH").
- `pairs`       : Dict of pairwise interaction parameters (e.g., "epsilon", "lambda" for square-well attractions).
- `species`     : Dict of user-defined species-level inputs (e.g., `input_densities`, `input_pSpecies`, `adjustable`).
- `system`      : Dict of global constants (e.g., "temperature", "bjerrum_length", etc.).
- `reactions`   : Vector of `ReactionEntry` objects parsed from the YAML input file.
- `bond_types`  : Vector of unique bonded bead-type index pairs, used in evaluating connectivity.
- `model`       : Vector of `Symbol`s indicating the free energy models requested (e.g., `[:ideal, :mfmt, :msa]`).
"""



struct PropertiesStruct
    monomers::Dict{Symbol, Any}              # e.g., "diameters" => [1.0, ...]
    pairs::Dict{Symbol, Any}                 # e.g., "square_well energys" => [1.0, ...]
    species::Dict{Symbol, Any}               # input densities, input pSpecies, adjustable
    system::Dict{Symbol, Float64}            # temperature, dielectric, etc.
    reactions::Vector{ReactionStruct}         # parsed from YAML
    bond_types::Vector{Tuple{Int, Int}}      # All unique bonded bead-type pairs
    fe_model :: Vector{Symbol}
end
