include("fourier_transforms/process_propagators.jl")
"""
    process_fourier(rho::BulkDensities, geometry::CartesianCoord) -> FourierCache

Initializes and returns a `FourierCache` containing preallocated Fourier-space arrays and FFTW plans
used for inhomogeneous DFT calculations.

# Arguments
- `rho::BulkDensities` — Density structure with bead and bond metadata.
- `geometry::CoordSystem` — Grid structure with system dimensions.

# Returns
- `FourierCache` — Struct holding transformed fields and FFT plans.
"""
function process_fourier(molsys::MolecularSystem, geometry::CartesianCoord)
    @unpack bin_width, NP, features = geometry
    @unpack mirrored, offset = features

    NP_star = compute_full_domain(NP, mirrored, offset)

    f_hat     = zeros(ComplexF64, NP_star)
    surf_hat  = zeros(ComplexF64, NP_star)

    # Monomer fields
    NB = length(molsys.properties.species[:monomers])
    dims_NB = (NP_star..., NB)
    rho_beads_hat    = zeros(ComplexF64, dims_NB)
    mu_ex_hat  = zeros(ComplexF64, dims_NB)

    # Bond fields
    MB = length(molsys.properties.bond_types)
    dims_MB = (NP_star..., 2, MB)
    rho_bonds_hat = zeros(ComplexF64, dims_MB)
    lng_hat       = zeros(ComplexF64, dims_MB)

    weight_bond_hat = process_propagator(molsys, geometry)

    planF = plan_fft!(f_hat, flags = FFTW.ESTIMATE)
    planB = plan_ifft!(f_hat, flags = FFTW.ESTIMATE)

    return FourierCache{ndims(surf_hat), ndims(mu_ex_hat), ndims(lng_hat),
                        typeof(planF), typeof(planB)}(
        rho_beads_hat, mu_ex_hat, rho_bonds_hat,
        weight_bond_hat, lng_hat,
        surf_hat, f_hat,
        planF, planB
    )
end