# === External Dependencies ===
using YAML
using UnPack

# === Core Data Structures ===
include("structures/MolecularSystem.jl")       # Top-level system container

# === Setup and Validation Logic ===
include("setup/io/check_yaml.jl")              # Input file format validation
include("setup/setup_molsys.jl")               # Constructs a MolecularSystem from a dataset

# === Public Constructor Function ===

"""
    molecular_system(input_file::String) -> MolecularSystem
    molecular_system(dataset::Dict) -> MolecularSystem

Construct a `MolecularSystem` from either a YAML file path or a pre-parsed dictionary.

If a file is provided, its structure is first validated using `check_yaml`.

# Arguments
- `input_file::String`: Path to a YAML input file.
- `dataset::Dict`: Pre-parsed YAML content.

# Returns
- `MolecularSystem`: Fully constructed system object with configurations, properties, and constants.
"""
function molecular_system(input_file::String)
    check_yaml(input_file)
    dataset = YAML.load_file(input_file)
    return setup_molsys(dataset)
end

function molecular_system(dataset::Dict)
    return setup_molsys(dataset)
end
