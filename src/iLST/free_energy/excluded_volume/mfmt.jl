"""
    struct MFMTFreeEnergy <: ExcessFreeEnergy

Stores the scalar and vector functional derivatives for the Modified Fundamental Measure Theory (MFMT)
hard-sphere free energy model in a bulk system.

# Fields
- `dPhi::Vector{Float64}`: Length-10 vector of functional derivatives ∂f/∂nₐ for scalar and vector weights.
- `nV1::Vector{Float64}`: 3-component vector weighted density nᵥ₁.
- `nV2::Vector{Float64}`: 3-component vector weighted density nᵥ₂.
"""
struct MFMTFreeEnergy <: ExcessFreeEnergy
    dPhi :: Vector{Float64}
    nV1  :: Vector{Float64}
    nV2  :: Vector{Float64}
end

"""
    construct_free_energy(::Type{MFMTFreeEnergy}, rho::BulkDensities) -> MFMTFreeEnergy

Initializes an `MFMTFreeEnergy` structure with zeroed fields appropriate for scalar and vector
weighted density derivatives used in MFMT calculations.

# Arguments
- `rho::BulkDensities`: Density structure (included for interface consistency).

# Returns
- `MFMTFreeEnergy`: Zero-initialized structure for MFMT model calculations.
"""
function construct_free_energy(::Type{MFMTFreeEnergy}, rho::BulkDensities)
    return MFMTFreeEnergy(
        zeros(10),  # dPhi
        zeros(3),   # nV1
        zeros(3)    # nV2
    )
end



"""
    compute_mfmt_scalars(ρ, DP) -> (n0, n1, n2, n3)

Computes scalar weighted densities for MFMT using the bulk bead density `ρ`
and diameter vector `DP`.
"""

function compute_mfmt_scalars(ρ, DP)
    n0 = sum(ρ)
    n1 = sum(@. @~ ρ * DP / 2)
    n2 = sum(@. @~ ρ * DP * DP * π)
    n3 = sum(@. @~ ρ * DP * DP * DP * π / 6)

    return n0, n1, n2, n3
end

"""
    free_energy(molsys :: MolecularSystem, bulk :: BulkState, fe_model::MFMTFreeEnergy) -> Float64

Computes the excess hard-sphere free energy using MFMT,
based on scalar weighted densities n₀, n₁, n₂, and n₃.
"""

function free_energy(molsys :: MolecularSystem, bulk :: BulkState, fe_model::MFMTFreeEnergy)
    rho  = bulk.rho.beads
    @unpack diameters = molsys.properties.monomers

    n0, n1, n2, n3 = compute_mfmt_scalars(rho, diameters)

    return free_energy_mfmt(n0, n1, n2, n3)
end


"""
    chemical_potential!(molsys :: MolecularSystem, bulk :: BulkState, fe_model::MFMTFreeEnergy)

Computes the MFMT excess chemical potential using scalar derivatives and
accumulates the result into `mu_ex` (bead-level chemical potential).
"""

function chemical_potential!(molsys :: MolecularSystem, bulk :: BulkState, fe_model :: MFMTFreeEnergy)
    
    dPhi = fe_model.dPhi
    
    rho  = bulk.rho.beads
    @unpack diameters = molsys.properties.monomers

    n0, n1, n2, n3 = compute_mfmt_scalars(rho, diameters)

    compute_mfmt_deriv!(dPhi, n0, n1, n2, n3)

    accumulate_mu_mfmt_bulk!(bulk.mu_ex, dPhi, diameters)
end


"""
    accumulate_mu_mfmt_bulk!(μ_ex, dΦ, σ)

Accumulates the MFMT excess chemical potential contribution into `μ_ex` for a bulk system,
based on scalar derivatives of the excess free energy functional.

# Arguments
- `μ_ex::Vector{Float64}`: Output vector of excess chemical potentials for each monomer/bead type. Modified in-place.
- `dΦ::AbstractVector{Float64}`: Vector of scalar functional derivatives 
  (∂fᵉˣ/∂n₀, ∂fᵉˣ/∂n₁, ∂fᵉˣ/∂n₂, ∂fᵉˣ/∂n₃).
- `σ::Vector{Float64}`: Vector of bead diameters, one per monomer type.

# Notes
This form applies only in the bulk where the weighted densities are computed from
scalar averages. The contributions are:

```math
μᵢ^{ex} += ∂f/∂n₀ + (∂f/∂n₁)(dᵢ/2) + (∂f/∂n₂)(π dᵢ²) + (∂f/∂n₃)(π dᵢ³ / 6)
"""

function accumulate_mu_mfmt_bulk!(μ_ex, dΦ :: AbstractVector, σ)
    @. μ_ex +=
            dΦ[1] +
            dΦ[2] * σ / 2 +
            dΦ[3] * π * σ^2 +
            dΦ[4] * π * σ^3 / 6
end

"""
    free_energy_mfmt(n0, n1, n2, n3, nV1, nV2) -> Float64

Computes the Modified Fundamental Measure Theory (MFMT) excess free energy density 
for a hard-sphere fluid using both scalar and vector weighted densities.

# Arguments
- `n0`, `n1`, `n2`, `n3::Float64`: Scalar weighted densities
- `nV1`, `nV2::SVector{3,Float64}` (or similar): Vector weighted densities

# Returns
- `f_hs::Float64`: The excess Helmholtz free energy density due to hard-sphere interactions.

# Notes
This full form includes directional correlations via dot products of `nV1` and `nV2`,
and is used in inhomogeneous systems such as density functional theory (DFT) grids.
"""


# inhomogeneous term
function free_energy_mfmt(n0, n1, n2, n3, nV1, nV2)
    one_minus_n3 = 1.0 - n3
    log_term     = log(one_minus_n3)
    inv_denom    = 1.0 / (36.0 * π * n3^2 * one_minus_n3^2)

    dot12 = dot(nV1, nV2)
    dot22 = dot(nV2, nV2)

    scalar_term = -n0 * log_term
    cross_term  = (n1 * n2 - dot12) / one_minus_n3
    cubic_term  = (n2^3 - 3.0 * n2 * dot22) * (n3 + one_minus_n3^2 * log_term) * inv_denom

    return scalar_term + cross_term + cubic_term
end


"""
    free_energy_mfmt(n0, n1, n2, n3) -> Float64

Computes the MFMT excess free energy density for a hard-sphere fluid using only scalar 
weighted densities. This simplified form omits vector contributions and is typically 
used for bulk systems.

# Arguments
- `n0`, `n1`, `n2`, `n3::Float64`: Scalar weighted densities

# Returns
- `f_hs::Float64`: Excess free energy density in the bulk (isotropic) approximation.

# Notes
Use this when directional correlations (e.g., from gradients or curvature) can be neglected.
"""

function free_energy_mfmt(n0, n1, n2, n3)
    one_minus_n3 = 1.0 - n3
    log_term     = log(one_minus_n3)
    inv_denom    = 1.0 / (36.0 * π * n3^2 * one_minus_n3^2)

    scalar_term = -n0 * log_term
    cross_term  = (n1 * n2) / one_minus_n3
    cubic_term  = n2^3 * (n3 + one_minus_n3^2 * log_term) * inv_denom

    return scalar_term + cross_term + cubic_term
end


"""
    compute_mfmt_deriv!(dΦ::AbstractVector, n0, n1, n2, n3, nV1, nV2)

Computes the scalar and vector derivatives of the MFMT (Modified Fundamental Measure Theory) 
hard-sphere excess free energy functional with respect to total weighted densities 
`n₀`, `n₁`, `n₂`, `n₃`, and vectorial contributions `nV₁`, `nV₂`.

# Arguments
- `dΦ::AbstractVector`: Output array of length 10. On return, contains:
    - Rows 1–4: ∂fᵉˣ/∂n₀, ∂fᵉˣ/∂n₁, ∂fᵉˣ/∂n₂, ∂fᵉˣ/∂n₃
    - Rows 5–7: ∂fᵉˣ/∂nV₁[1–3]
    - Rows 8–10: ∂fᵉˣ/∂nV₂[1–3]
- `n0`, `n1`, `n2`, `n3::Float64`: Scalar weighted densities
- `nV1`, `nV2::SVector{3,Float64}`: Vector weighted densities

# Notes
This version is used for systems where the MFMT functional depends on total weighted densities,
typically in bulk or single-component inhomogeneous systems.
"""


function compute_mfmt_deriv!(dΦ::AbstractVector, n0, n1, n2, n3, nV1, nV2)

    one_minus_n3      = 1.0 - n3
    log_1_minus_n3    = log(one_minus_n3)
    inv_1_minus_n3    = 1.0 / one_minus_n3
    inv_n3            = 1.0 / n3
    inv_36pi          = 1.0 / (36.0 * π)
    
    dot12 = dot(nV1, nV2)
    dot22 = dot(nV2, nV2)

    # Scalar terms
    dΦ[1] = -log_1_minus_n3

    dΦ[2] = n2 * inv_1_minus_n3

    dΦ[3] = n1 * inv_1_minus_n3 +
                   (n2^2 - dot22) * 3.0 * inv_36pi * inv_n3 *
                   (log_1_minus_n3 * inv_n3 + inv_1_minus_n3^2)

    dΦ[4] = n0 * inv_1_minus_n3 +
                   (n1 * n2 - dot12) * inv_1_minus_n3^2 -
                   (log_1_minus_n3 * inv_n3^3 * 2.0 * inv_36pi +
                    inv_n3^2 * inv_1_minus_n3 * inv_36pi +
                    (1.0 - 3.0 * n3) * inv_n3^2 * inv_36pi * inv_1_minus_n3^3) *
                   (n2^3 - 3.0 * n2 * dot22)

    # Vector terms
    for v in 1:3
        dΦ[v + 4] = -nV2[v] * inv_1_minus_n3
        dΦ[v + 7] = -nV1[v] * inv_1_minus_n3 -
                           (log_1_minus_n3 * inv_n3 + inv_1_minus_n3^2) *
                           n2 * nV2[v] * inv_n3 * 6.0 * inv_36pi
    end
end

"""
    compute_mfmt_deriv!(dΦ, n0, n1, n2, n3)

Computes the scalar derivatives of the Modified Fundamental Measure Theory (MFMT) 
hard-sphere excess free energy with respect to the four scalar weighted densities:
n₀, n₁, n₂, and n₃.

These correspond to:
- `n0`: total number density
- `n1`: surface area-like contribution
- `n2`: curvature-related contribution
- `n3`: volume packing fraction

Assumes a homogeneous (bulk) system, so vector dot products `dot12` and `dot22` are zero.

# Arguments
- `dΦ::Vector{Float64}`: Output array (length 4) to store ∂Fᵉˣ/∂nᵢ for i = 0 to 3
- `n0`, `n1`, `n2`, `n3::Float64`: Scalar weighted densities of the system

# Modifies
- `dΦ` is updated in-place with the four partial derivatives.

# Reference
Rosenfeld, Y. (1989). Free-energy model for the inhomogeneous hard-sphere fluid mixture and density-functional theory of freezing. *Phys. Rev. Lett.*, 63(9), 980.

"""

function compute_mfmt_deriv!(dΦ :: AbstractVector, n0, n1, n2, n3)

    one_minus_n3      = 1.0 - n3
    log_1_minus_n3    = log(one_minus_n3)
    inv_1_minus_n3    = 1.0 / one_minus_n3
    inv_n3            = 1.0 / n3
    inv_36pi          = 1.0 / (36.0 * π)
    
    # For bulk systems (vector contributions are zero)
    dot12 = 0.0
    dot22 = 0.0
    
    dΦ[1] = -log_1_minus_n3
    
    dΦ[2] = n2 * inv_1_minus_n3
    
    dΦ[3] = n1 * inv_1_minus_n3 +
                   (n2^2 - dot22) * 3.0 * inv_36pi * inv_n3 *
                   (log_1_minus_n3 * inv_n3 + inv_1_minus_n3^2)
    
    dΦ[4] = n0 * inv_1_minus_n3 +
                   (n1 * n2 - dot12) * inv_1_minus_n3^2 -
                   (log_1_minus_n3 * inv_n3^3 * 2.0 * inv_36pi +
                    inv_n3^2 * inv_1_minus_n3 * inv_36pi +
                    (1.0 - 3.0 * n3) * inv_n3^2 * inv_36pi * inv_1_minus_n3^3) *
                   (n2^3 - 3.0 * n2 * dot22)
    
end