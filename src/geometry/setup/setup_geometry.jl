include("coordinate_system/create_coordsys.jl")

"""
    setup_geometry(dataset::Dict) -> CartesianGrid{D}

Parses the `geometry` section from the input dataset and constructs the appropriate grid
object (currently supports Cartesian only). This function is typically called during system
initialization.

Returns:
- `CartesianGrid{D}`: the spatial discretization grid.
"""

function setup_geometry(dataset::Dict)
    return create_coordinate_system(dataset) 
end
