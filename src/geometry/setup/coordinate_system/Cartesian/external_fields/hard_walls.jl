"""
    cartesian_features(features::Dict, ::Val{:hard_walls})

Processes the `:hard_walls` external field during `GeometrySystem` setup.

This function:
- Sets `periodic[dim] = false` in each spatial direction where hard walls are specified.
- Initializes `position` vectors for each specified dimension, marking hard wall locations as `[0.0, box_length]`.

This function does **not** evaluate numerical external fields — that is handled later by the DFT solver (e.g., via `eval_External`).

Expected input format (inside `dataset["geometry"]["external_field"]`):
    external_field:
      hard_walls:
        dims: [1, 3]  # 1-based indices for dimensions with hard walls

Modifies:
- `features[:periodic]`                             → sets to `false` in specified directions
- `features[:external_field][:hard_walls][:position]` → sets default wall positions
"""

struct ExtHardWalls <: AbstractExternalField 
    positions :: Vector{Vector{Float64}}
end

function cartesian_features(features::CartesianFeatures, data_external_field, ::Val{:hard_walls})
   @unpack dimensions, periodic, mirrored = features

    description = data_external_field[:hard_walls]

    if !haskey(description, "dims")
        error("Missing 'dims' key in external_field[:hard_walls]. Please specify affected directions.")
    end

    dims = description["dims"]
    ndims = length(periodic)

    # Initialize position array: one entry per direction
    description["position"] = [Float64[] for _ in 1:ndims]

    for dim in dims
        @assert 1 ≤ dim ≤ ndims "Invalid dimension index: $dim (expected between 1 and $ndims)"

        # Turn off periodicity
        periodic[dim] = false

        # Define default wall positions at 0.0 and box edge
        description["position"][dim] = [0.0, dimensions[dim]]
    end

    return ExtHardWalls(description["position"])
end

"""
    update_features(features::Dict, ::Val{:hard_walls})

Stub for `:hard_walls` external field. Currently no updates are needed after mirroring.
"""
function update_features(features::CartesianFeatures, data_external_field, ::Val{:hard_walls})
    # Nothing to update
end
