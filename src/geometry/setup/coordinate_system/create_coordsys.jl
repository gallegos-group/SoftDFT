include("Cartesian/create_cartesian.jl")

"""
    create_coordinate_system(dataset::Dict) -> CartesianGrid{D}

Constructs a coordinate system grid from the `geometry` section of the dataset.

Supported values for `coordinate_system`:
- "Cartesian" → returns `CartesianGrid{D}` (D = 1, 2, or 3)

Throws an error if the coordinate system is missing or unsupported.

Example:

    dataset = Dict(
        "geometry" => Dict(
            "coordinate_system" => "Cartesian",
            "dimensions" => [10.0, 10.0, 20.0],
            "bin_width" => [0.5, 0.5, 1.0]
        )
    )
    grid = create_coordinate_system(dataset)
"""
function create_coordinate_system(dataset::Dict)
    geom = get(dataset, "geometry", nothing)
    isnothing(geom) && error("Missing required section: 'geometry'.")

    coord_symbol = get(geom, "coordinate_system", nothing)
    isnothing(coord_symbol) && error("Missing required key: 'coordinate_system' under 'geometry'.")

    system = Symbol(coord_symbol)

    if system == :Cartesian
        return create_cartesian(dataset)
    elseif system == :Spherical
        error("Spherical coordinate system not yet supported.")
    elseif system == :Cylindrical
        error("Cylindrical coordinate system not yet supported.")
    else
        error("Unknown coordinate system: $system. Supported: :Cartesian.")
    end
end
