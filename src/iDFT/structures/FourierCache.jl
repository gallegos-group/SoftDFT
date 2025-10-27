"""
    FourierCache{M, N, O, F, B}

Stores Fourier-space representations of grid-based fields and FFT plans for efficient
evaluation of convolution and field operations in inhomogeneous DFT.

This structure contains forward and backward FFT plans, preallocated arrays for storing
transformed densities, excess potentials, and interaction kernels, and temporary storage 
used in electrostatics or surface-related computations.

# Type Parameters
- `M` — Dimensionality of surface-related FFT arrays.
- `N` — Dimensionality of monomer-related FFT arrays.
- `O` — Dimensionality of bond-related FFT arrays.
- `F`, `B` — Types of the forward and backward FFT plans, subtypes of `AbstractFFTs.Plan`.

# Fields
- `rho_beads_hat::Array{ComplexF64, N}` — Fourier transform of monomer densities.
- `mu_ex_hat::Array{ComplexF64, N}` — Fourier transform of excess chemical potentials.

- `rho_bonds_hat::Array{ComplexF64, O}` — Fourier transform of bond densities.
- `weight_bonds_hat::Array{Float64, N}` — Fourier-space weights for bonds.
- `lng_hat::Array{ComplexF64, O}` — Fourier transform of logarithmic bond weights.

- `surf_hat::Array{ComplexF64, M}` — Surface-related field or kernel in Fourier space.
- `f_hat::Array{ComplexF64, M}` — Temporary field for general-purpose FFT operations.

- `plan_forward::F` — FFT plan for forward transforms.
- `plan_backward::B` — FFT plan for inverse transforms.
"""
struct FourierCache{M, N, O, F<:AbstractFFTs.Plan, B<:AbstractFFTs.Plan}
    rho_beads_hat     :: Array{ComplexF64, N}
    mu_ex_hat         :: Array{ComplexF64, N}

    mu_sim_hat         :: Array{ComplexF64, N}

    rho_bonds_hat     :: Array{ComplexF64, O}
    weight_bonds_hat  :: Array{Float64, N}
    lng_hat           :: Array{ComplexF64, O}

    surf_hat          :: Array{ComplexF64, M}
    f_hat             :: Array{ComplexF64, M}

    plan_forward      :: F
    plan_backward     :: B
end
