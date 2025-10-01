"""
    parse_chemistry(dataset::Dict) 
        -> (species_dict::Dict, reactions::Vector{ReactionEntry}, 
            used_monomers::Dict{String, Int}, state_families::Dict{String, Vector{String}})

Parse and organize all chemistry-related input fields from the YAML dataset.

This includes:
- Extracting raw species data and constructing the `species_dict` that holds names, densities, pSpecies values, and adjustability flags.
- Parsing chemical reactions and formatting them into `ReactionEntry` structs.
- Identifying all monomers used in sequences or reactions.
- Building and validating state families for protonation or tautomerization, ensuring consistency and proper disjoint partitioning.

Returns:
- `species_dict`: Dictionary of species fields needed for further processing
- `reaction_entries`: List of parsed reaction stoichiometries and pK values
- `used_monomers`: Mapping from monomer name to index used in sequences
- `state_families`: Mapping from base-state monomer to list of all states in that family
"""


include("chemistry/parse_species.jl")
include("chemistry/parse_reactions.jl")
include("chemistry/used_monomers.jl")
include("chemistry/state_families.jl")

function parse_chemistry(dataset::Dict, constants :: ConstantsStruct)
    # Parse reaction entries from input YAML
    reaction_entries = parse_reactions(dataset)

    # Extract raw species data
    species_data = dataset["species"]

    species_entries = parse_species(dataset)

    used_monomers  = build_used_monomers(species_data, reaction_entries)

    state_families = build_state_families(species_entries, reaction_entries)

    # Build species dict fields
    species_names     = [s["sequence"] for s in species_data]
    
    input_densities = Float64[]
    for s in species_data
        if haskey(s, "concentration")
            conc = s["concentration"]
            ρ = concentration_to_density(conc, constants)
            push!(input_densities, ρ)
        else
            ρ = get(s, "density", 0.0)
            push!(input_densities, ρ)
        end
    end
    
    input_pSpecies    = [get(s, "pSpecies", NaN) for s in species_data]
    adjustable_flags  = [get(s, "adjustable", false) for s in species_data]
    compensator_flags = [get(s, "adjustable", false) for s in species_data]
    species_charge    = [0.0 for s in species_data]
    role = [get(s, "role", :context) for s in species_data]

    species_dict = Dict(
        "species"         => species_names,
        "input_densities" => input_densities,
        "input_pSpecies"  => input_pSpecies,
        "adjustable"      => adjustable_flags,
        "compensator"     => compensator_flags,
        "species_charge"  => species_charge,
        "monomers"        => sort(collect(keys(used_monomers))),
        "state_families"  => state_families,
        "role" =>  role
    )

    return species_dict, reaction_entries, used_monomers, state_families
end
