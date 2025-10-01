include("process_cartesian.jl")

"""
    create_cartesian(dataset::Dict) -> CartesianGrid{D}

Constructs a `CartesianGrid{D}` object from the provided `dataset`, which must contain
a `geometry` section with the keys `dimensions` and `bin_width`.

# Required Fields in `dataset["geometry"]`
- `dimensions::Vector{Float64}`: The physical size of the domain in each coordinate direction.
- `bin_width::Vector{Float64}`: The width of each bin in the corresponding direction.

# Optional Fields
- `external_field::Dict`: Defines spatially varying external fields
- `features::Dict{String, Any}`: Optional user-provided flags; these are merged with derived geometry features.

# Returns
- A `CartesianGrid{D}` instance with dimensions, bin widths, and feature metadata.

# Throws
- An error if dimensions and bin widths are not specified, mismatched, or unsupported.
"""
function create_cartesian(dataset::Dict)
    geom = dataset["geometry"]

    if !haskey(geom, "dimensions") || !haskey(geom, "bin_width")
        error("Error: Please specify both `dimensions` and `bin_width` under `geometry`.")
    end

    bins_vec = geom["bin_width"]
    D = length(bins_vec)

    if D != length(geom["dimensions"])
        error("Error: Mismatch in number of `dimensions` and `bin_width` entries.")
    elseif D ∉ (1, 2, 3)
        error("Error: Only 1D, 2D, or 3D Cartesian grids are supported.")
    end

    bin_width  = Tuple{Vararg{Float64, D}}(bins_vec)

    # Use process_cartesian to modify dimensions + build features
    dimensions, features = process_cartesian(dataset, bin_width)

    # Convert to NTuple form
    NP         = ntuple(i -> Int(cld(dimensions[i], bin_width[i])) + 1, D)

    return CartesianGrid{D}(dimensions, bin_width, NP, features)
end