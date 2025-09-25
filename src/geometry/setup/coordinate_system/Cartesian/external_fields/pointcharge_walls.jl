"""
    cartesian_features(features::Dict, ::Val{:pointcharge_walls})

Processes the `:pointcharge_walls` external field during `GeometrySystem` setup.

This function:
- Marks specified directions as non-periodic.
- Initializes wall positions as `[0.0, box_length]` for each affected dimension.
- Stores left/right point charges in `features[:external_field][:pointcharge_walls][:charges]`.
- Updates `features[:total_charge]` accordingly.
- Sets `mirrored[dim] = true` if the left/right point charges differ.

This function does **not** apply the external field — that is handled by the DFT solver via `eval_External`.

Expected input format:
    external_field:
      pointcharge_walls:
        dims: [1, 2]
        charge: [[-0.05, -0.05], [0.05, 0.07]]  # qL and qR per dimension

Modifies:
- `features[:periodic]`
- `features[:mirrored]`
- `features[:external_field][:pointcharge_walls][:position]`
- `features[:external_field][:pointcharge_walls][:charges]`
- `features[:total_charge]`
"""
function cartesian_features(features::Dict{Symbol, Any}, ::Val{:pointcharge_walls})
    description = features[:external_field][:pointcharge_walls]

    dims = get(description, "dims", nothing)
    isnothing(dims) && error("Please specify the wall dimensions under 'dims'.")

    dimensions = features[:dimensions]
    periodic   = features[:periodic]
    mirrored   = features[:mirrored]

    charge = get(description, "charge", nothing)
    isnothing(charge) && error("Please provide wall charges under 'charge'.")

    # Normalize charge format
    if charge isa Vector{Float64}
        charge = [charge]
    elseif !(charge isa Vector{<:AbstractVector{<:Real}})
        error("Charge must be a Vector of Vectors, e.g., [[qL, qR], [qL, qR], ...].")
    end

    if length(charge) != length(dims)
        error("Mismatch between number of 'dims' and 'charge' entries.")
    end

    D = length(dimensions)
    description["charge"] = charge
    description["position"] = [Float64[] for _ in 1:D]
    description["charges"]  = [Float64[] for _ in 1:D]

    for (i, dim) in enumerate(dims)
        @assert 1 ≤ dim ≤ D "Invalid dimension index: $dim"

        periodic[dim] = false
        qL, qR = charge[i]

        description["position"][dim] = [0.0, dimensions[dim]]
        description["charges"][dim]  = [qL, qR]

        features[:total_charge][dim] += qL + qR

        if qL != qR
            mirrored[dim] = true
        end
    end
end

function update_features(features::Dict{Symbol, Any}, ::Val{:pointcharge_walls})
    # description = features[:external_field][:pointcharge_walls]
    # dims = description["dims"]
    # mirrored = features[:mirrored]

    # for (i, dim) in enumerate(dims)
    #     if mirrored[dim]
    #         qL, qR = description["charge"][i]
    #         features[:total_charge][dim] += qL + qR
    #     end
    # end
end # I dont think this needs to be here