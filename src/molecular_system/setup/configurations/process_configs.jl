"""
    process_configurations(
        dataset::Dict,
        used_monomers::Dict{String, Int},
        state_families::Dict{String, Vector{String}},
        properties::PropertiesStruct,
    ) -> Vector{ConfigurationStruct}

Builds the configuration objects for each species listed in the input dataset. This routine constructs:

- `sequence`: a list of monomer indices used to represent the polymer.
- `states`: a vector indicating the number of protonation or chemical states per bead (default 1).
- `state_family`: for each bead, lists all indices in its allowed state family.
- `topology`: the hierarchical structural representation of the species (linear, branched, etc.).
- `species_model`: the physical model (e.g., freely jointed chain) as determined by species options.

All `ConfigurationStruct` objects are collected and returned as a vector. This forms the core of the molecular system description used for simulation.
"""


include("species_model.jl")
include("parse_topology.jl")




function process_configurations(
    dataset::AbstractDict,
    used_monomers::AbstractDict,
    state_families::AbstractDict,
)
    species_data = dataset["species"]
    configurations = ConfigurationStruct[]

    # Build state index lookup: monomer => state index
    state_index = Dict{String, Int}()
    for family in values(state_families)
        sorted_family = sort(family)
        for (idx, monomer) in enumerate(sorted_family)
            state_index[monomer] = idx
        end
    end

    for species in species_data

        sequence_str = species["sequence"]
        topology = parse_topology(sequence_str)
        n_beads = length(topology.monomers)

        species_model = determine_species_model(species, n_beads)

        sequence = Vector{Int}(undef, n_beads)
        states = ones(Int, n_beads)
        state_family = Vector{Vector{Int}}(undef, n_beads)

        for i in 1:n_beads
            monomer_name = string(topology.monomers[i])
            monomer_index = get(used_monomers, monomer_name, nothing)
        
            if monomer_index === nothing
                error("Monomer $monomer_name not found in used_monomers.")
            end
        
            sequence[i] = monomer_index
        
            if haskey(state_families, monomer_name)
                # Convert the state family from names to indices
                name_family = state_families[monomer_name]
                index_family = [used_monomers[n] for n in name_family]
                states[i] = length(index_family)
                state_family[i] = index_family
            else
                states[i] = 1
                state_family[i] = [monomer_index]
            end
        end

        push!(configurations, ConfigurationStruct(species_model, sequence, states, topology, state_family))
    end
    
    return configurations
end

