"""
    cartesian_features(features::Dict, ::Val{:sw_walls})

Processes the `:sw_walls` (square-well wall potential) external field during `GeometrySystem` setup.

This function:
- Sets `periodic[dim] = false` for all specified dimensions.
- Initializes default wall positions as `[0.0, box_length]` in each specified direction.
- Stores the wall positions in `features[:external_field][:sw_walls][:position]`.

Note: This only sets up metadata; the actual potential is evaluated later by the DFT solver.

Expected input format:
    external_field:
      sw_walls:
        dims: [1, 2]

Modifies:
- `features[:periodic]`
- `features[:external_field][:sw_walls][:position]`
"""
function cartesian_features(features::Dict{Symbol, Any}, ::Val{:sw_walls})
    dimensions = features[:dimensions]
    periodic   = features[:periodic]
    mirrored   = features[:mirrored]
    field_spec = features[:external_field][:sw_walls]

    dims = get(field_spec, "dims", nothing)
    isnothing(dims) && error("Missing required key: 'dims' under 'sw_walls'. Please specify as 'dims: [i, j, ...]'.")

    n_dims = length(dimensions)
    field_spec["position"] = [Float64[] for _ in 1:n_dims]

    for dim in dims
        if dim < 1 || dim > n_dims
            error("Invalid dimension $dim in 'sw_walls.dims'. Must be between 1 and $n_dims.")
        end
        periodic[dim] = false
        field_spec["position"][dim] = [0.0, dimensions[dim]]
    end
end

"""
    update_features(features::Dict, ::Val{:sw_walls})

Stub for `:sw_walls` external field. No updates are required for mirrored domains.
"""
function update_features(features::Dict{Symbol, Any}, ::Val{:sw_walls})
    # No update needed for now
end