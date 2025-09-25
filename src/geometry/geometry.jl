# === External Dependencies ===
using YAML
using UnPack

# === Coordinate System Abstract Type ===
abstract type CoordSystem end

# === Core Data Structures ===
include("structures/Cartesian.jl")

# === Setup Logic ===
include("setup/setup_geometry.jl")

# === Public Constructor Function ===

"""
    geometry_system(input_file::String) -> GeometrySystem
    geometry_system(dataset::Dict)       -> GeometrySystem

Construct a `GeometrySystem` from either a YAML file path or a pre-parsed dictionary.

- If a file path is provided, the file is loaded using `YAML.load_file`.
- If a dictionary is provided, it is used directly.

Returns a `GeometrySystem` (currently, this is just a coordinate grid like `CartesianGrid{D}`).
"""
function geometry_system(input_file::String)
    dataset = YAML.load_file(input_file)
    return setup_geometry(dataset)
end

function geometry_system(dataset::Dict)
    return setup_geometry(dataset)
end
