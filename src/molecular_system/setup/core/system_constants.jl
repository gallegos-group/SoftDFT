"""
    system_constants(dataset::Dict) -> (ConstantsStruct, Dict{String, Float64})

Extract physical constants and global system-level parameters from the YAML dataset.

This function performs:
- Parsing of user-defined or default reduced units (e.g., molecular length scale `l_sc`)
- Extraction of temperature and dielectric constant
- Computation of the Bjerrum length using fundamental constants
- Assembly of a user-facing dictionary of system properties

# Arguments
- `dataset::Dict`: Parsed YAML input dictionary

# Returns
- `ConstantsStruct`: Contains physical constants and reduced unit scaling
- `Dict{String, Float64}`: Dictionary of system-level properties, including:
    - `"temperature"`
    - `"dielectric_constant"'
    - `"bjerrum_length"`
"""

# --- Dependency includes ---
include("constants/process_constants.jl")         # Parse reduced units like l_sc
include("conversions/conversions.jl")         # Parse reduced units like l_sc
include("system/get_temperature.jl")           # Extract temperature from dataset
include("system/get_dielectric_constant.jl")   # Extract dielectric constant
include("system/get_Bjerrumlength.jl")         # Compute Bjerrum length
include("system/process_system.jl")                # Assemble final system property dictionary

function system_constants(dataset::Dict)
    constants = process_constants(dataset)
    system    = process_system(dataset, constants)
    return constants, system
end
