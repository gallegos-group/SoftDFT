"""
    setup_molsys(input_file::String) -> MolecularSystem
    setup_molsys(dataset::Dict) -> MolecularSystem

Top-level function to initialize a `MolecularSystem` from a YAML input file or parsed dataset.

This function performs all steps required to construct a self-consistent molecular system,
including:

- Validating the input file structure.
- Parsing physical constants and system-level properties (e.g., temperature, dielectric constant).
- Extracting species, reactions, and state transition networks.
- Computing monomer-level properties and reaction free energy shifts (Δμ).
- Generating polymer configurations with sequence data, topologies, and bonding constraints.

# Arguments
- `input_file::String`: Path to the YAML input file.
- `dataset::Dict`: (Advanced use) A pre-parsed YAML dictionary.

# Returns
A `MolecularSystem` struct containing:
- `configurations::Vector{ConfigurationStruct}`: All user-defined polymers or molecules.
- `properties::PropertiesStruct`: Encodes monomer/pair interactions, species data, and reactions.
- `constants::ConstantsStruct`: Physical constants including reduced unit scaling.

# Example
```julia
system = setup_molsys("input.yaml")
This function serves as the main entry point for downstream DFT or simulation workflows.
"""

# === Core Physical Constants & System-Wide Properties ===
include("core/system_constants.jl")           # get_temperature, get_dielectric_constant, etc.
include("core/parse_free_energy.jl")      # parse free energy model list

# === Species & Chemical Reactions ===
include("species/parse_chemistry.jl")     # reactions, used monomers, state families
include("species/parse_properties.jl")    # monomer properties, delta_muH

# === Polymer Configurations ===
include("configurations/process_configs.jl")   # parse sequence + states
include("configurations/process_bonds.jl")     # extract bond types and topologies

# ---------------------------
# === Main Setup Function ===
# ---------------------------

function setup_molsys(dataset::Dict)

    constants, system_properties = system_constants(dataset)

    fe_model = parse_free_energy(dataset["model"])

    species_dict, reaction_entries, used_monomers, state_families = 
                parse_chemistry(dataset, constants)

    configurations = process_configurations(dataset, used_monomers, state_families)

    bond_types = process_bonds(configurations)

    properties = PropertiesStruct(  Dict{Symbol, Float64}(),
                                    Dict{Symbol, Any}(),
                                    SpeciesProperties(symbolize_keys(species_dict)),
                                    symbolize_keys(system_properties),
                                    reaction_entries,
                                    bond_types,
                                    fe_model)

    parse_species_properties(dataset, used_monomers, state_families,
                                constants, configurations, properties)

    bonding_constraint(configurations, properties)

    return MolecularSystem(configurations, properties, constants)
end

function symbolize_keys(d::Dict)
    Dict(Symbol(k) => v isa Dict ? symbolize_keys(v) : v for (k, v) in d)
end

function stringify_keys(d::Dict)
    Dict(string(k) => v isa Dict ? stringify_keys(v) : v for (k, v) in d)
end