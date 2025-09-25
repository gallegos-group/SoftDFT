"""
    parse_species_properties(
        dataset::Dict,
        used_monomers::Dict{String, Int},
        state_families::Dict{String, Vector{String}},
        constants::ConstantsStruct,
        properties::PropertiesStruct,
    ) -> Nothing

Parse and assign all physical and thermodynamic properties associated with monomers and their interactions.

This routine:
- Loads physical descriptors for each monomer (e.g., diameter, valence) using `process_monomers`.
- Computes the ideal work (ΔμH) for protonation/deprotonation transitions based on reactions and species data via `assign_delta_muH!`.
- Identifies pairwise interaction parameters (e.g., square-well interaction strengths and ranges) through `process_pairs`.

All parsed data are stored in-place within the `properties.monomers` and `properties.pairs` dictionaries.
"""


include("properties/process_monomers.jl")
include("properties/assign_delta_muH.jl")
include("properties/process_pairs.jl")
include("properties/assign_species_charge.jl")

function parse_species_properties(dataset, used_monomers, state_families, constants, configurations, properties)

    process_monomers(used_monomers, dataset, properties)

    assign_delta_muH!(properties, used_monomers, state_families, constants)

    process_pairs(used_monomers, dataset, properties)

    assign_species_charge!(configurations, properties)
end