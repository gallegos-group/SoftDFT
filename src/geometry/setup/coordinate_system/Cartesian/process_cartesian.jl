include("cartesian_features.jl")

"""
    process_cartesian(dataset::Dict) -> NTuple{D, Float64}, Dict{Symbol, Any}

Processes the `geometry` section of the input dataset to prepare geometry dimensions
and field-related features, including offsets, periodicity, and external field effects.

Returns the updated dimensions (after offset and mirroring) and the feature dictionary.
"""



function process_cartesian(dataset::Dict, dims_vec, bins_vec)
    geometry_data = dataset["geometry"]
    D = length(dims_vec)

    # Base geometry metadata
    dimensions    = copy(dims_vec)
    periodic      = fill(true, D)
    mirrored      = fill(true, D)
    offset        = zeros(Int, D)
    total_charge  = zeros(Float64, D)

    # --- Apply offset for non-periodic directions ---
    offset_value = 1.0
    for i in 1:D
        if !periodic[i]
            offset[i] = round(Int, offset_value / bins_vec[i])
        end
    end

    # Construct typed features struct
    features = CartesianFeatures(dimensions, periodic, mirrored, offset, total_charge)

    # Optional external field definitions
    data_external_field = haskey(geometry_data, "external_field") ?
       Dict(Symbol(k) => v for (k, v) in geometry_data["external_field"]) : Dict()

    # Apply feature-specific logic (sets periodicity, positions, charges, etc.)
    external_field = cartesian_features(features, data_external_field)

    @unpack offset, periodic, mirrored, total_charge = features
    return dimensions, external_field, offset, periodic, mirrored, total_charge, features
end
