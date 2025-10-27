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

# auto-detect version
function create_cartesian(dataset::Dict)
    D = length(dataset["geometry"]["dimensions"])
    return D == 1 ? create_cartesian(dataset, Val(1)) :
           D == 2 ? create_cartesian(dataset, Val(2)) :
           D == 3 ? create_cartesian(dataset, Val(3)) :
           error("Unsupported dimensionality $D")
end

function create_cartesian(dataset::Dict, ::Val{D}) where {D}
    geom = dataset["geometry"]

    if !haskey(geom, "dimensions") || !haskey(geom, "bin_width")
        error("Error: Please specify both `dimensions` and `bin_width` under `geometry`.")
    end

    bins_raw = geom["bin_width"]
    dims_raw = geom["dimensions"]

    bins_vec = Float64[]
    dims_vec = Float64[]

    for x in bins_raw
        push!(bins_vec, Float64(x))
    end

    for x in dims_raw
        push!(dims_vec, Float64(x))
    end

    dims = length(bins_vec)

    if dims != length(dims_vec)
        error("Error: Mismatch in number of `dimensions` and `bin_width` entries.")
    elseif dims ∉ (1, 2, 3)
        error("Error: Only 1D, 2D, or 3D Cartesian grids are supported.")
    end

    # Use process_cartesian to modify dimensions + build features
    dimensions, external_field, offset, periodic, mirrored, total_charge, features = process_cartesian(dataset, dims_vec, bins_vec)

    # Convert to NTuple form
    NP = ntuple(i -> Int(cld(dimensions[i], bins_vec[i])) + 1, D)

    return CartesianGrid(Tuple(dimensions), Tuple(bins_vec), NP, periodic, mirrored, offset, total_charge, external_field, features)
end