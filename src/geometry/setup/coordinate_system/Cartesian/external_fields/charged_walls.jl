"""
    cartesian_features(features::Dict, ::Val{:charged_walls})

Processes the `:charged_walls` external field during GeometrySystem setup.

This function:
- Sets `periodic[dim] = false` in each specified direction.
- Initializes default wall positions as `[0.0, box_length]` for each charged wall.
- Stores left/right surface charge densities in `features[:external_field][:charged_walls][:charges]`.
- Updates `features[:total_charge][dim]` by summing wall charges.
- Sets `mirrored[dim] = true` if the left/right charges are unequal.

This function does **not** modify the actual external field array — that is done later by the DFT solver (e.g., via `eval_External`).

Expected input format (inside `dataset["geometry"]["external_field"]`):
    external_field:
      charged_walls:
        dims: [1, 3]
        charge: [[-0.01, -0.01], [0.02, 0.02]]  # charge[i] = [left, right] for dim i

Modifies:
- `features[:periodic]`
- `features[:mirrored]`
- `features[:external_field][:charged_walls][:position]`
- `features[:external_field][:charged_walls][:charges]`
- `features[:total_charge]`
"""

struct ExtChargedWalls <: AbstractExternalField 
    positions :: Vector{Vector{Float64}}
    charges :: Vector{Vector{Float64}}
end

function cartesian_features(features::CartesianFeatures, data_external_field, ::Val{:charged_walls})
    @unpack dimensions, periodic, mirrored, total_charge = features

    description = data_external_field[:charged_walls]

    # Validate and extract dimensions to apply wall charges
    dims = get(description, "dims", nothing)
    isnothing(dims) && error("Please specify the dimensions using 'dims: [i, j, ...]'.")

    D = length(dimensions)
    @assert all(1 .≤ dims .≤ D) "Invalid dimension index in dims: $dims"

    # Normalize charge input to Vector{Vector{Float64}}
    charge = get(description, "charge", nothing)
    isnothing(charge) && error("Please specify wall charges under 'charge'.")

    if charge isa Vector{Float64}
        charge = [charge]  # promote to vector of vectors
    elseif !(charge isa Vector{<:AbstractVector{<:Real}})
        error("Charge must be a vector of vectors, e.g., [[qL₁, qR₁], [qL₂, qR₂], ...].")
    end

    if length(charge) != length(dims)
        error("Mismatch: number of charge entries does not match dims.\nDims = $dims, Charge = $charge")
    end

    # Initialize and populate geometry-specific data
    description["position"] = [Float64[] for _ in 1:D]
    description["charges"]  = [Float64[] for _ in 1:D]

    for (i, dim) in enumerate(dims)
        periodic[dim] = false
        qL, qR = charge[i]

        description["position"][dim] = [0.0, dimensions[dim]]
        description["charges"][dim]  = [qL, qR]

        total_charge[dim] += qL + qR

        if qL != qR
            mirrored[dim] = true
        end
    end

    # Normalize and store charge data back
    description["charge"] = charge

    return ExtChargedWalls(description["position"], description["charges"])
end

"""
    update_features(features::Dict, ::Val{:charged_walls})

Updates `features[:total_charge]` for mirrored directions by adding the left and right
charges a second time, accounting for domain doubling due to symmetry.

Only applies to directions where `mirrored[dim] == true`.
"""
function update_features(features::CartesianFeatures, data_external_field, ::Val{:charged_walls})
    # description = features[:external_field][:charged_walls]
    # dims = description["dims"]
    # mirrored = features[:mirrored]

    # for (i, dim) in enumerate(dims)
    #     if mirrored[dim]
    #         qL, qR = description["charge"][i]
    #         features[:total_charge][dim] += qL + qR
    #     end
    # end
end # I dont think this needs to be here
