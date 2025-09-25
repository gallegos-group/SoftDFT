"""
    cartesian_features(features::Dict, ::Val{:input_external})

Processes the `:input_external` field specification for reading external potential
data from an external file (e.g., CSV, text grid). Does not modify geometry features,
but validates structure and initializes `:position` for compatibility.

Expected input format:
    external_field:
      input_external:
        dims: [1, 2, 3]  # directions the field applies to
        file: "path/to/potential.dat"  # (optional) external file to be read later

Notes:
- This function only prepares geometry metadata and validates input.
- Actual field loading and application should be handled later (e.g., in `eval_External`).
- You may choose to enable setting `periodic` or `mirrored` in the future.

Modifies:
- `features[:external_field][:input_external][:position]` → initialized as empty per dimension
"""
function cartesian_features(features::Dict{Symbol, Any}, ::Val{:input_external})
    description = features[:external_field][:input_external]

    dims = get(description, "dims", nothing)
    isnothing(dims) && error("Error: please specify the dimensions using 'dims: [i, j, ...]'.")

    n_dims   = length(features[:dimensions])
    periodic = features[:periodic]
    mirrored = features[:mirrored]

    # Initialize positions array for compatibility (even if unused here)
    description["position"] = [Float64[] for _ in 1:n_dims]

    # Optional: allow this field to affect periodicity or mirroring if needed
    for dim in dims
        @assert 1 ≤ dim ≤ n_dims "Invalid dimension index: $dim"
        # Optional logic:
        # periodic[dim] = false
        # mirrored[dim] = false
    end
end

"""
    update_features(features::Dict, ::Val{:input_external})

Stub for `:input_external`. No updates are currently required after mirroring.
"""
function update_features(features::Dict{Symbol, Any}, ::Val{:input_external})
    # Nothing to be done
end
