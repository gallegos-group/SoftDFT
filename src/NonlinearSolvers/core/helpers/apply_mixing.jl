"""
    apply_mixing(previous, current, amix) -> new_object

Blend `previous` and `current` using linear mixing:
    (1 - amix) * previous + amix * current

### Behavior
- Returns a new object (non-mutating).
- Use this for scalars or tuples of scalars/arrays.
- Recursively applies mixing to tuple elements.

### Arguments
- `previous`, `current`: Scalars or tuples of compatible structure.
- `amix::Float64`: Mixing parameter (typically ∈ [0, 1]).

### Returns
- A new value or tuple of values.

Use `apply_mixing!` if in-place modification is desired.
"""
function apply_mixing(previous::Number, current::Number, amix::Float64)
    return (1.0 - amix) * previous + amix * current
end

function apply_mixing(previous::Tuple, current::Tuple, amix::Float64)
    return ntuple(i -> apply_mixing(previous[i], current[i], amix), length(previous))
end

"""
    apply_mixing!(previous, current, amix)

Overwrite `previous` with a linear blend of `previous` and `current`:
    previous .= (1 - amix) * previous + amix * current

### Behavior
- Modifies `previous` in-place.
- Supports arrays and nested tuples of arrays.
- Errors if any element in a tuple is a scalar (immutable).

Use `apply_mixing` for scalar or read-only contexts.
"""
@inline function apply_mixing!(previous::AbstractArray,
                               current::AbstractArray,
                               amix::Float64)
    coeff1 = 1.0 - amix
    coeff2 = amix
    @inbounds @simd for i in eachindex(previous)
        previous[i] = coeff1 * previous[i] + coeff2 * current[i]
    end
    return previous
end

@inline function apply_mixing!(previous::AbstractVector{<:AbstractArray},
                               current::AbstractVector{<:AbstractArray},
                               amix::Float64)
    @inbounds for j in eachindex(previous)
        apply_mixing!(previous[j], current[j], amix)
    end
    return previous
end



# function apply_mixing!(previous::Tuple, current::Tuple, amix::Float64)
#     for i in eachindex(previous)
#         if previous[i] isa Number
#             error("Cannot modify tuple field $i (a Number) in-place. Use `apply_mixing` instead.")
#         end
#         apply_mixing!(previous[i], current[i], amix)
#     end
# end

@inline function apply_mixing!(p::Tuple{}, c::Tuple{}, amix::Float64)
    return p
end

@inline function apply_mixing!(p::Tuple, c::Tuple, amix::Float64)
    @inbounds apply_mixing!(p[1], c[1], amix)
    @inbounds apply_mixing!(Base.tail(p), Base.tail(c), amix)
    return p
end

# Optional fallback to catch unsupported types early
@noinline function apply_mixing!(previous, current, amix::Float64)
    error("apply_mixing! not defined for type: $(typeof(previous))")
end
