"""
    FixedSpecies

Container for representing a species with one or more *fixed* segments in inhomogeneous DFT calculations.
This is typically used when certain molecules or configurations are immobilized or externally imposed
(e.g., grafted chains, adsorbed species, or surface-bound segments).

# Fields

- `density::Vector{Float64}` — Total fixed density per species (often scalar-valued per species).
- `segments::Vector{Bool}` — Boolean vector indicating which segments are fixed (`true` if fixed).
- `configuration::Vector{Array{Float64, 2}}` — Predefined spatial density profiles for each fixed segment.
- `coordinates::Vector{Vector{Tuple}}` — Coordinate lists for each fixed segment (e.g., grid positions).
- `lambda::Vector{Float64}` — Scaling factors for each fixed species, used in normalization or reweighting.

This structure is used during initialization and evolution of systems where some parts of the molecule
remain spatially constrained or serve as external inputs to the system.
"""

struct FixedSpecies{D}
    density       :: Vector{Float64}                 # Total fixed density per species
    segments      :: Vector{Bool}                    # Marks which segments are fixed
    coordinates   :: Vector{NTuple{D, Int}}          # Coordinates for each segment (per config)
    lambda        :: Vector{Float64}                 # Weight or scaling factor for fixed species
end

"""
    FixedSpecies(num_segments::Int) -> FixedSpecies

Create a `FixedSpecies` object with all fields initialized for a species with `num_segments` segments.
"""
function FixedSpecies(num_segments::Int, ndims :: Int)
    density       = [0.0]
    segments      = fill(false, num_segments)
    coordinates   = [ntuple(_ -> 0, ndims) for _ in 1:num_segments]
    lambda        = [1.0]

    return FixedSpecies(density, segments, coordinates, lambda)
end
