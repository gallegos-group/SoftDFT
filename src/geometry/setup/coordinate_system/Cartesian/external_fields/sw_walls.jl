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

struct ExtSquareWellWalls <: AbstractExternalField
    names :: Vector{String}
    positions :: Vector{Vector{Float64}}
    energys :: Vector{Vector{Vector{Float64}}}
    lambdas :: Vector{Vector{Vector{Float64}}}
end

# function cartesian_features(features::CartesianFeatures, data_external_field, ::Val{:sw_walls})
#     @unpack dimensions, periodic, mirrored = features

#     field_spec = data_external_field[:sw_walls]

#     dims = get(field_spec, "dims", nothing)
#     isnothing(dims) && error("Missing required key: 'dims' under 'sw_walls'. Please specify as 'dims: [i, j, ...]'.")

#     n_dims = length(dimensions)
#     field_spec["position"] = [Float64[] for _ in 1:n_dims]

#     for dim in dims
#         if dim < 1 || dim > n_dims
#             error("Invalid dimension $dim in 'sw_walls.dims'. Must be between 1 and $n_dims.")
#         end
#         periodic[dim] = false
#         field_spec["position"][dim] = [0.0, dimensions[dim]]
#     end

#     # --- Parse monomer-specific energy/lambda data ---
#     names = Vector{String}()
#     energys = Vector{Vector{Float64}}()
#     lambdas = Vector{Vector{Float64}}()

#     if haskey(field_spec, "monomers")
#         monomers = field_spec["monomers"]
#         for (monomer, data) in monomers
#             push!(names, monomer)
#             energy_data = get(data, "energy", nothing)
#             lambda_data = get(data, "lambda", nothing)
#             if isnothing(energy_data) || isnothing(lambda_data)
#                 error("Each monomer entry under 'sw_walls.monomers' must define both 'energy' and 'lambda'.")
#             end

#             # Convert [[a,b]] → Vector{Float64}
#             push!(energys, [Float64(x) for x in energy_data[1]])
#             push!(lambdas, [Float64(x) for x in lambda_data[1]])
#         end
#     else
#         error("Missing required section: 'monomers' under 'sw_walls'.")
#     end

#     return ExtSquareWellWalls(names, field_spec["position"], energys, lambdas)
# end

function cartesian_features(features::CartesianFeatures, data_external_field, ::Val{:sw_walls})
    @unpack dimensions, periodic, mirrored = features

    field_spec = data_external_field[:sw_walls]

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

    # --- Parse monomer-specific energy/lambda data ---
    names = Vector{String}()
    monomer_energy = Dict{String, Vector{Vector{Float64}}}()
    monomer_lambda = Dict{String, Vector{Vector{Float64}}}()

    if haskey(field_spec, "monomers")
        monomers = field_spec["monomers"]
        for (monomer, data) in monomers
            push!(names, monomer)

            energy_data = get(data, "energy", nothing)
            lambda_data = get(data, "lambda", nothing)

            if isnothing(energy_data) || isnothing(lambda_data)
                error("Each monomer entry under 'sw_walls.monomers' must define both 'energy' and 'lambda'.")
            end

            # Normalize to per-dimension form
            # e.g. [[-3.0, -0.5]] → replicate for all n_dims if only one provided
            if length(energy_data) == 1
                energy_data = fill(energy_data[1], n_dims)
                lambda_data = fill(lambda_data[1], n_dims)
            elseif length(energy_data) != n_dims
                error("Energy/lambda for monomer $monomer must either be length 1 or length $n_dims.")
            end

            monomer_energy[monomer] = [Float64.(vec) for vec in energy_data]
            monomer_lambda[monomer] = [Float64.(vec) for vec in lambda_data]
        end
    else
        error("Missing required section: 'monomers' under 'sw_walls'.")
    end

    # --- Repack into dimension-major order ---
    energys  = [ [ monomer_energy[name][dim] for name in names ] for dim in 1:n_dims ]
    lambdas  = [ [ monomer_lambda[name][dim] for name in names ] for dim in 1:n_dims ]

    return ExtSquareWellWalls(names, field_spec["position"], energys, lambdas)
end

"""
    update_features(features::Dict, ::Val{:sw_walls})

Stub for `:sw_walls` external field. No updates are required for mirrored domains.
"""
function update_features(features::CartesianFeatures, data_external_field, ::Val{:sw_walls})
    # No update needed for now
end