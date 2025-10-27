include("DensityFields.jl")
include("ExcessFields.jl")
include("FourierCache.jl")
include("FixedSpecies.jl")

"""
    SpatialFields

Container for all spatially-resolved field variables required in inhomogeneous DFT calculations.

This structure bundles density fields, excess potentials, Fourier transform caches, and
information on fixed species. It serves as the primary interface for evaluating external 
fields, updating densities, and computing convolution-based interactions.

# Fields
- `rho_K::DensityFields` — Bead and bond densities, including simulated components.
- `excess::ExcessFields` — External, excess chemical, electrostatic, and integration fields.
- `fourier::FourierCache` — Cached forward/backward FFT plans and transformed arrays.
- `fixed::Vector{FixedSpecies}` — Per-configuration data for grafted or fixed segments.
"""

struct SpatialFields{D, E, F, FS}
    rho_K      :: D             # DensityFields
    excess     :: E             # ExcessFields
    fourier    :: F             # FourierCache
    fixed      :: FS            # Vector{FixedSpecies}
end
