"""
    abstract type CartesianCoord <: CoordSystem

Abstract supertype for all Cartesian coordinate representations.

---

    struct CartesianGrid{D} <: CartesianCoord

Represents a regular Cartesian grid in D dimensions.

# Fields
- `dimensions :: NTuple{D, Float64}`: Physical size of the domain in each direction.
- `bin_width  :: NTuple{D, Float64}`: Width of each bin (grid spacing).
- `NP         :: NTuple{D, Int}`: Number of bins in each direction.
- `features   :: Dict{Symbol, Any}`: Metadata including boundary conditions, wall specs, etc.

---

# Convenient Type Aliases
- `CartesianZ`   → 1D grid
- `CartesianYZ`  → 2D grid
- `CartesianXYZ` → 3D grid
"""

# === Cartesian Coordinate System Abstract Type ===
abstract type CartesianCoord <: CoordSystem end

# === Cartesian Grid Structure ===
struct CartesianGrid{D} <: CartesianCoord
    dimensions :: NTuple{D, Float64}      # Domain size in each direction
    bin_width  :: NTuple{D, Float64}      # Bin width in each direction
    NP         :: NTuple{D, Int}          # Number of bins in each direction
    features   :: Dict{Symbol, Any}       # Optional tags (e.g., symmetry, wall info)
end

# === Convenient Type Aliases ===
const CartesianZ   = CartesianGrid{1}
const CartesianYZ  = CartesianGrid{2}
const CartesianXYZ = CartesianGrid{3}