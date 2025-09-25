function apply_mirroring!(field_K::AbstractArray, geometry)
    NP       = geometry.NP                   # true dimensions
    offset   = geometry.features[:offset]    # ::NTuple
    mirrored = geometry.features[:mirrored]  # ::Vector{Bool}
    NDIM     = length(NP)

    spatial_iter  = CartesianIndices(NP)
    trailing_iter = CartesianIndices(size(field_K)[NDIM+1:end])

    for K in spatial_iter
        K_star = to_star_index(K, offset, mirrored)

        for J in trailing_iter
            val = @inbounds field_K[K_star, J]

            # assign into all mirrored partners (including itself)
            for mirrored_K in mirrored_index_iterator(K_star, NP, offset, mirrored)
                # if K != mirrored_K println(Tuple(K)[1],"     ",Tuple(mirrored_K)[1]) end
                @inbounds field_K[mirrored_K, J] = val
            end
        end
    end
    # Something is going wrong here
    return field_K
end


"""
    to_star_index(K::CartesianIndex, offsets::NTuple{N, Int}, mirrored::Vector{Bool})

Convert true–domain index `K` into calculation–domain index `K_star`
by adding offsets. For mirrored dimensions, this gives the
location in the primary half.
"""
@inline function to_star_index(K::CartesianIndex{N},
                               offsets::Vector{Int},
                               mirrored::Vector{Bool}) where {N}
    @assert length(offsets) == N
    @assert length(mirrored) == N
    return CartesianIndex(ntuple(i -> K[i] + offsets[i], N))
end


"""
    from_star_index(K_star::CartesianIndex, NP::NTuple{N,Int},
                    offsets::NTuple{N, Int}, mirrored::Vector{Bool})

Map a calc–domain index back into the true–domain index.
For mirrored dims, folds values from the mirrored half
onto the base interval.
"""
@inline function from_star_index(K_star::CartesianIndex{N},
                                 NP::NTuple{N, Int},
                                 offsets::Vector{Int},
                                 mirrored::Vector{Bool}) where {N}
    @assert length(offsets) == N
    @assert length(mirrored) == N

    return CartesianIndex(ntuple(i -> begin
        k = K_star[i] - offsets[i]
        if mirrored[i]
            # fold into [1, NP[i]]
            L = NP[i]
            k = if k < 1
                1
            elseif k > L
                2L - k + 1
            else
                k
            end
        end
        k
    end, N))
end
