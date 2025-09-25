"""
    DensityFields{M, N}

Container for storing density-related fields used in inhomogeneous DFT calculations.

# Fields
- `beads::Array{Float64, M}` — Density of individual beads (monomers) across a discretized spatial grid.
- `bonds::Array{Float64, N}` — Density of bonded pairs (e.g., for evaluating chain connectivity) across the grid.
- `simulated::Array{Float64, M}` — Simulated or accumulated bead densities from particle configurations.
- segment densities

# Type Parameters
- `M` — Dimensionality of bead-related arrays (e.g., `(Nx, Ny, Nz, Nb)`).
- `N` — Dimensionality of bond-related arrays (e.g., `(Nx, Ny, Nz, 2, Nbond)`).

Typically constructed during setup and updated throughout the DFT solver to track spatial and species-resolved densities.
"""

struct DensityFields{M, N}
    beads     :: Array{Float64, M}
    bonds     :: Array{Float64, N}
    simulated :: Array{Float64, M}
    segments  :: Vector{Array{Float64, N}}
end