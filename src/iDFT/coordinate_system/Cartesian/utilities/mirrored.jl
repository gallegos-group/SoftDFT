"""
    mirrored_index_iterator(K_star::CartesianIndex,
                            NP::NTuple{N,Int},
                            offset::NTuple{N,Int},
                            mirrored::Vector{Bool})

Yield all mirrored indices in storage space (including K_star itself).
For each mirrored dimension, reflection is done about the plane
centered between `lo = offset[d] + 1` and `hi = offset[d] + NP[d]`.
"""

function mirrored_index_iterator(K_star::CartesianIndex{N},
                                 NP::NTuple{N,Int},
                                 offset::Vector{Int},
                                 mirrored::Vector{Bool}) where {N}
    idx = Tuple(K_star)
    mirror_states = Iterators.product(ntuple(_ -> (false, true), N)...)
    NP_star = compute_full_domain(NP, mirrored, offset)
    return Iterators.map(mirror_state -> begin
        reflected = ntuple(d -> mirrored[d] && mirror_state[d] ? NP_star[d] - idx[d] + 1 : idx[d], N)
        CartesianIndex(reflected)
    end, mirror_states)
end

# function mirrored_index_iterator(K_star::CartesianIndex{N},
#                                  NP::NTuple{N,Int},
#                                  offset::Vector{Int},
#                                  mirrored::Vector{Bool}) where {N}
#     idx = Tuple(K_star)

#     mirror_states = Iterators.product(ntuple(_ -> (false, true), N)...)

#     return Iterators.map(mirror_state -> begin
#         reflected = ntuple(d -> begin
#             if mirrored[d] && mirror_state[d]
#                 lo = offset[d] + 1
#                 hi = offset[d] + NP[d]
#                 lo + hi - idx[d]
#             else
#                 idx[d]
#             end
#         end, N)
#         CartesianIndex(reflected)
#     end, mirror_states)
# end

"""
    compute_full_domain(NP::NTuple{N, Int}, mirrored::Vector{Bool}, offsets::Vector{Int})

Compute the effective domain size `NP_star` from the true domain `NP`.

- `NP`: true domain sizes (without offsets or mirroring).
- `mirrored`: flags for each dimension indicating whether it is mirrored.
- `offsets`: number of ghost/offset points to add on each side of the dimension.

The result includes both offsets and mirroring.
"""
function compute_full_domain(NP::NTuple{N, Int}, mirrored::Vector{Bool}, offsets::Vector{Int}) where N
    @assert length(mirrored) == N "Length of 'mirrored' must match length of NP"
    @assert length(offsets) == N "Length of 'offsets' must match length of NP"
    return ntuple(i -> (NP[i] + 2*offsets[i]) * (mirrored[i] ? 2 : 1), N)
end
