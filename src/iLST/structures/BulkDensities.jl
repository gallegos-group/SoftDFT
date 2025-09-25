"""
    BulkDensities

Container for all density fields used in the bulk thermodynamic calculation.

Fields:
- `species`  :: Vector of total densities for each species.
- `segments` :: Vector of matrices. `segments[i]` is a matrix of size (n_states, n_segments)
                giving state-resolved densities for each segment of species `i`.
- `pairs`    :: Vector of 4D arrays. `pairs[i][j_state, j_seg, i_state, i_seg]` gives the
                pairwise density between two segment-state combinations within species `i`.
- `beads`    :: Vector of bead (monomer) densities aggregated over all species and states.
- `bonds`    :: 2×n matrix of bond densities. `bonds[1, b]` and `bonds[2, b]` correspond to
                left and right contributions to bond type `b`, respectively.
"""

struct BulkDensities
    species  :: Vector{Float64}                          # One value per species
    segments :: Vector{Matrix{Float64}}                  # segments[i][state, segment]
    pairs    :: Vector{Array{Float64, 4}}                # pairs[i][state_j, segment_j, state_i, segment_i]
    beads    :: Vector{Float64}                          # One value per unique monomer/bead type
    bonds    :: Matrix{Float64}                          # One value per unique bond
end