include("SpatialFields.jl")
include("../functionals/Abstract_Types/AbstractFunctional.jl")

"""
    IsingDFT

Encapsulates all components required to solve an inhomogeneous density functional theory (DFT)
problem for weak polyelectrolyte systems using the Ising DFT framework.

This structure is constructed from a bulk reference state (`IsingLST`), a coordinate grid, 
a list of excess free energy functionals, and spatially resolved field variables. It serves 
as the core data object passed through setup, evaluation, and analysis routines for iDFT.

# Fields
- `bulk_system::IsingLST` — Precomputed bulk state (molecular system and thermodynamic reference).
- `geometry::CoordSystem` — Spatial grid and coordinate representation (Cartesian, spherical, etc.).
- `fields::SpatialFields` — Spatial density, excess, FFT, and fixed segment data.
- `functionals::Vector{AbstractFunctional}` — Collection of free energy contributions (e.g., ideal, MF, SW).
- `numerics::Dict{String, Any}` — Numerical parameters for convergence, tolerances, and solver behavior.
"""

struct IsingDFT
    bulk_system :: IsingLST
    geometry :: CoordSystem
    fields :: SpatialFields
    functionals :: Vector{AbstractFunctional}
    numerics :: Dict{String, Any}
end
