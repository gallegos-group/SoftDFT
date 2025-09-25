include("cartesian_features.jl")

"""
    process_cartesian(dataset::Dict) -> NTuple{D, Float64}, Dict{Symbol, Any}

Processes the `geometry` section of the input dataset to prepare geometry dimensions
and field-related features, including offsets, periodicity, and external field effects.

Returns the updated dimensions (after offset and mirroring) and the feature dictionary.
"""
function process_cartesian(dataset::Dict, bin_width)
    geometry_data = dataset["geometry"]
    dims_vec = geometry_data["dimensions"]
    D = length(dims_vec)

    # Copy base dimensions
    dimensions = copy(dims_vec)
    features = Dict{Symbol, Any}()

    # Initialize base geometry metadata
    features[:dimensions] = copy(dimensions)
    features[:periodic]   = fill(true, D)
    features[:mirrored]   = fill(false, D)
    features[:offset]     = zeros(Int, D)
    features[:total_charge] = zeros(Float64, D)

    # Optional external field definitions
    features[:external_field] = haskey(geometry_data, "external_field") ?
        Dict(Symbol(k) => v for (k, v) in geometry_data["external_field"]) : Dict()


    # Apply feature-specific logic (sets periodicity, positions, charges, etc.)
    cartesian_features(features)

    # --- Apply offset for non-periodic directions ---
    offset_value = 1.0  # must be > convolution width
    for i in 1:D
        if !features[:periodic][i]
            features[:offset][i] = round(Int, offset_value / bin_width[i])
        end
    end

    return Tuple(dimensions), features
end
