"""
    apply_flat_solution!(solution_variables::Tuple, X_k::Vector{Float64}) -> Nothing

Applies a flat vector of values `X_k` to a structured collection of `solution_variables`,
updating them in-place. This is the inverse of `flatten_solution!` and is typically used to
map optimization or solver updates back into structured simulation variables.

### Arguments
- `solution_variables::Tuple`: A nested structure composed of mutable containers, such as:
  - `AbstractArray` (e.g. `Array{Float64, N}`)
  - `AbstractVector{<:AbstractArray}` (e.g. vector of 1D or 2D arrays)
  - Nested `Tuple`s of the above
  - Optionally `ContinuationVariable`, assuming `set_value!` is defined
- `X_k::Vector{Float64}`: A flat vector whose values will overwrite the contents of `solution_variables`.
  Must be the same length as the result of `flatten_solution!(..., solution_variables)`.

### Behavior
- The structure and shapes of `solution_variables` are preserved.
- Values are read from `X_k` in order and written into each container using `CartesianIndices`, supporting arbitrary array dimensions.
- Dispatch-based design allows clean extension to new types.

### Returns
- `Nothing`. The `solution_variables` are updated in-place.

### Throws
- An error if an unsupported type is encountered in `solution_variables`.

### Example
```julia
solution = (rand(2, 2), [rand(3), rand(3)])
X_k = zeros(10)
flatten_solution!(X_k, solution)  # Save original state
X_k .= 42.0                       # Overwrite with dummy data
apply_flat_solution!(solution, X_k)  # Push flat vector back into structured solution
Notes

    solution_variables must consist only of mutable types.

    Scalars are not supported.

    Assumes consistency between the number of entries in X_k and the structure of solution_variables.
"""

# Entry point: starts recursive application from index 1
function apply_flat_solution!(solution_variables::Tuple, X_k::Vector{Float64})
    _ = apply_flat_solution!(solution_variables, X_k, 1)
    return nothing
end

# Case 1: Tuple: recurse through each element
function apply_flat_solution!(vars::Tuple, X_k::Vector{Float64}, idx::Int)
    for v in vars
        idx = apply_flat_solution!(v, X_k, idx)
    end
    return idx
end

# Case 2: Vector of arrays
function apply_flat_solution!(var::AbstractVector{<:AbstractArray}, X_k::Vector{Float64}, idx::Int)
    for u in eachindex(var), K in CartesianIndices(var[u])
        var[u][K] = X_k[idx]
        idx += 1
    end
    return idx
end

# Case 3: Single array
function apply_flat_solution!(var::AbstractArray, X_k::Vector{Float64}, idx::Int)
    @inbounds for K in CartesianIndices(var)
        var[K] = X_k[idx]
        idx += 1
    end
    return idx
end

# Optional: ContinuationVariable (assumes set_value! is defined)
# if @isdefined ContinuationVariable
#     function apply_flat_solution!(var::ContinuationVariable, X_k::Vector{Float64}, idx::Int)
#         set_value!(var, X_k[idx])
#         return idx + 1
#     end
# end

# Fallback for unsupported types
function apply_flat_solution!(var, X_k::Vector{Float64}, idx::Int)
    error("Unsupported type in solution_variables: $(typeof(var))")
end
