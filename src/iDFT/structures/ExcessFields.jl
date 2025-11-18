"""
    ExcessFields{L, M, N}

Holds the excess field contributions used in inhomogeneous DFT calculations, including
external potentials, excess chemical potentials, bond weights, and electrostatic fields.

# Fields
- `mu_ex_K::Array{Float64, M}` — Excess chemical potential field for each bead across the spatial grid.
- `lng_K::Array{Float64, N}` — Logarithmic bonding field used in chain connectivity contributions.
- `Ext::Array{Float64, M}` — External field contributions (e.g., hard walls, soft walls, point charges).
- `trapez::Array{Float64, M}` — Trapezoidal rule weights for spatial integration.
- `Psi::Array{Float64, L}` — Electrostatic potential field over the grid.
- `PsiC::Vector{Float64}` — Mean correction term for the electrostatic potential (e.g., to ensure neutrality).

# Type Parameters
- `L` — Dimensionality of the scalar potential `Psi` (e.g., `(Nx, Ny, Nz)`).
- `M` — Dimensionality of bead-level fields (e.g., `(Nx, Ny, Nz, Nb)`).
- `N` — Dimensionality of bond-level fields (e.g., `(Nx, Ny, Nz, 2, Nbond)`).

Used to accumulate and apply non-ideal interactions and external driving fields in the DFT solver.
"""

struct ExcessFields{L, M, N}
    lambda_K :: Array{Float64, M}   # one-body field
    mu_ex_K :: Array{Float64, M}    # excess chemical potential on grid
    lng_K   :: Array{Float64, N}    # logarithmic bond weight
    Ext     :: Array{Float64, M}    # soft/hard wall or field-like term
    trapez  :: Array{Float64, M}    # trapezoidal integration weights
    Psi     :: Array{Float64, L}    # electrostatic potential field
    PsiC    :: Vector{Float64}      # mean correction to electrostatics
end
