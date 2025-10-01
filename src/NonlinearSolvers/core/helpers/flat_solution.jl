"""
    flatten_solution!(X_k::Vector{Float64}, solution_variables::Tuple)

Flattens a structured collection of solution variables into a 1D vector of `Float64` values,
writing them into `X_k` in-place. This function supports deeply nested and heterogeneous structures
via multiple dispatch and is typically used to prepare a structured solution for linear solvers,
mixing schemes, or continuation methods.

### Arguments
- `X_k::Vector{Float64}`: A preallocated flat vector to store the contents of `solution_variables`.
  Must be large enough to hold all flattened values.
- `solution_variables::Tuple`: A nested structure composed of arrays, or vectors of arrays.
  Scalars are not supported. Valid elements include:
  - `AbstractArray` (e.g., `Array{Float64, N}`)
  - `AbstractVector{<:AbstractArray}` (e.g., array-of-arrays)
  - Nested `Tuple`s of the above

### Method Internals
This function delegates to multiple dispatch methods of the form:
```julia
flatten_solution!(X_k::Vector{Float64}, var, idx::Int)
"""

# Entry point: starts recursive flattening from index 1
function flatten_solution!(X_k::Vector{Float64}, solution_variables :: Tuple)
    final_idx = flatten_solution!(X_k, solution_variables, 1)
    return nothing
end

# Tuple: recurse through each element
function flatten_solution!(X_k::Vector{Float64}, vars::Tuple, idx::Int)
    for v in vars
        idx = flatten_solution!(X_k, v, idx)
    end
    return idx
end

# Vector of arrays (e.g., multiple field contributions)
function flatten_solution!(X_k::Vector{Float64}, var::AbstractVector{<:AbstractArray}, idx::Int)
    for u in eachindex(var), K in CartesianIndices(var[u])
        X_k[idx] = var[u][K]
        idx += 1
    end
    return idx
end

# A single array
function flatten_solution!(X_k::Vector{Float64}, var::AbstractArray, idx::Int)
    @inbounds for K in CartesianIndices(var)
        X_k[idx] = var[K]
        idx += 1
    end
    return idx
end

# Fallback: error if unknown type
function flatten_solution!(X_k::Vector{Float64}, var, idx::Int)
    error("Unsupported type in solution_variables: $(typeof(var))")
end