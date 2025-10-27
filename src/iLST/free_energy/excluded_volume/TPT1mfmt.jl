"""
    struct TPT1MFMTFreeEnergy <: ExcessFreeEnergy

Holds state and intermediate arrays required to compute the MFMT + TPT1 free energy 
and its functional derivatives in bulk systems.

# Fields
- `dPhi::Matrix{Float64}`: (10 × N) matrix of functional derivatives ∂f/∂nₐ for each bead type.
- `lng::Vector{Float64}`: Logarithm of hard-sphere contact values y_HS for each bond type.
- `n0i`, `n1i`, `n2i`, `n3i`::Vector{Float64}: Per-bead scalar weighted densities.
- `nV1`, `nV2`::Vector{Float64}: Global vectorial weighted densities (3 components each).
- `wt_bonds::Vector{Float64}`: Optional work vector for bond-specific intermediates.
"""
struct TPT1MFMTFreeEnergy <: ExcessFreeEnergy
    dPhi     :: Matrix{Float64}
    lng      :: Vector{Float64}
    n0i      :: Vector{Float64}
    n1i      :: Vector{Float64}
    n2i      :: Vector{Float64}
    n3i      :: Vector{Float64}
    nV1      :: Vector{Float64}
    nV2      :: Vector{Float64}
    wt_bonds :: Matrix{Float64}
end


"""
    construct_free_energy(::Type{TPT1MFMTFreeEnergy}, rho::BulkDensities) -> TPT1MFMTFreeEnergy

Initializes a zeroed `TPT1MFMTFreeEnergy` structure with appropriate dimensions
based on the provided density structure `rho`.

# Arguments
- `rho::BulkDensities`: Contains `.beads` and `.bonds` arrays for determining field sizes.

# Returns
- `TPT1MFMTFreeEnergy`: Struct with allocated and zero-initialized fields.
"""
function construct_free_energy(::Type{TPT1MFMTFreeEnergy}, rho::BulkDensities)
    nbeads = length(rho.beads)
    nbonds = size(rho.bonds, 2)

    return TPT1MFMTFreeEnergy(
        zeros(10, nbeads),  # dPhi
        zeros(nbonds),      # lng
        zeros(nbeads),      # n0i
        zeros(nbeads),      # n1i
        zeros(nbeads),      # n2i
        zeros(nbeads),      # n3i
        zeros(3),           # nV1
        zeros(3),           # nV2
        zeros(2, nbonds)       # wt_bonds
    )
end



    
"""
    compute_mfmt_scalars_i!(n0i, n1i, n2i, n3i, ρ, σ)

Computes the scalar weighted densities for each species in the MFMT functional.

# Arguments
- `n0i, n1i, n2i, n3i::Vector{Float64}`: Output vectors to store weighted densities.
- `ρ::Vector{Float64}`: Bead density for each species.
- `σ::Vector{Float64}`: Diameter of each species.

# Notes
Each output is computed as:
- `n₀ᵢ = ρᵢ`
- `n₁ᵢ = ρᵢ * σᵢ / 2`
- `n₂ᵢ = ρᵢ * π * σᵢ²`
- `n₃ᵢ = ρᵢ * π * σᵢ³ / 6`
"""
function compute_mfmt_scalars_i!(n0i, n1i, n2i, n3i, ρ, σ)
    @. n0i = ρ
    @. n1i = ρ * σ / 2
    @. n2i = ρ * π * σ^2
    @. n3i = ρ * π * σ^3 / 6
end


"""
    free_energy(molsys, bulk, fe_model) -> Float64

Computes the excess free energy for a bulk system using the Modified Fundamental Measure Theory (MFMT)
with chain connectivity corrections via the Thermodynamic Perturbation Theory (TPT1).

# Arguments
- `molsys::MolecularSystem`: Molecular system data (e.g. diameters, bond types).
- `bulk::BulkState`: Bulk densities of beads and bonds.
- `fe_model::TPT1MFMTFreeEnergy`: Temporary storage for weighted densities.

# Returns
- Scalar excess Helmholtz free energy density (in k_BT units).
"""
function free_energy(
    molsys  :: MolecularSystem,
    bulk    :: BulkState,
    fe_model:: TPT1MFMTFreeEnergy
)

    rho        = bulk.rho.beads
    rho_bonds  = bulk.rho.bonds

    @unpack bond_types = molsys.properties
    @unpack diameters = molsys.properties.monomers

    @unpack n0i, n1i, n2i, n3i = fe_model

    compute_mfmt_scalars_i!(n0i, n1i, n2i, n3i, rho, diameters)

    return free_energy_TPT1mfmt(rho_bonds, bond_types, n0i, n1i, n2i, n3i, diameters)
end

"""
    chemical_potential(molsys, bulk, fe_model)

Computes the MFMT + TPT1 excess chemical potential for a bulk system.

# Arguments
- `molsys::MolecularSystem`: Molecular system containing monomer diameters and bond types.
- `bulk::BulkState`: Bulk state containing bead and bond densities, and where μᵉˣ is stored.
- `fe_model::TPT1MFMTFreeEnergy`: Storage for weighted densities, log(y_HS), and derivatives.

# Modifies
- `bulk.mu_ex`: Updated in-place with MFMT + TPT1 excess chemical potential.
- `bulk.lng`: Updated in-place with log(y_HS) from TPT1 correction.
"""

function chemical_potential!(
    molsys  :: MolecularSystem,
    bulk    :: BulkState,
    fe_model:: TPT1MFMTFreeEnergy
)

    rho = bulk.rho.beads
    rho_bonds = bulk.rho.bonds
    
    @unpack bond_types = molsys.properties
    @unpack diameters = molsys.properties.monomers

    @unpack n0i, n1i, n2i, n3i, lng, dPhi = fe_model

    compute_mfmt_scalars_i!(n0i, n1i, n2i, n3i, rho, diameters)

    compute_TPT1mfmt_deriv!(dPhi, lng, rho_bonds, bond_types, n0i, n1i, n2i, n3i, diameters)

    accumulate_mu_mfmt_bulk!(bulk.mu_ex, dPhi, diameters)

    @. bulk.lng += lng
end


"""
    accumulate_mu_mfmt_bulk!(mu_ex, dΦ, σ)

Accumulates the MFMT excess chemical potential contribution into `mu_ex` for a bulk system,
based on species-resolved scalar derivatives of the excess free energy functional.

# Arguments
- `mu_ex::Vector{Float64}`: Output vector of excess chemical potentials for each monomer/bead type (modified in-place).
- `dΦ::AbstractMatrix{Float64}`: Matrix of scalar functional derivatives, size (4, N), 
  where each column corresponds to one bead type:
    - Row 1: ∂fᵉˣ/∂n₀
    - Row 2: ∂fᵉˣ/∂n₁
    - Row 3: ∂fᵉˣ/∂n₂
    - Row 4: ∂fᵉˣ/∂n₃
- `σ::Vector{Float64}`: Vector of bead diameters, one per bead type.

# Notes
This is used in bulk systems where only scalar weighted densities are defined. The excess chemical potential for each bead `i` is given by:

```math
μ_i^{ex} = ∂f/∂n₀ + (∂f/∂n₁)(σᵢ/2) + (∂f/∂n₂)(π σᵢ²) + (∂f/∂n₃)(π σᵢ³ / 6)
"""

function accumulate_mu_mfmt_bulk!(mu_ex, dΦ :: AbstractMatrix, σ)
    @views @. mu_ex +=
                        dΦ[1, :] +
                        dΦ[2, :] * σ / 2 +
                        dΦ[3, :] * π * σ^2 +
                        dΦ[4, :] * π * σ^3 / 6
end

"""
    free_energy_TPT1mfmt(rho_bonds, bond_types, n0i, n1i, n2i, n3i, nV1, nV2, σ)

Computes the excess free energy density for a multi-component hard-sphere fluid using 
Modified Fundamental Measure Theory (MFMT) and TPT1 chain connectivity corrections,
including vectorial (inhomogeneous) correlations.

# Arguments
- `rho_bonds::Vector{Float64}`: Density of each bonded pair.
- `bond_types::Vector{Tuple{Int, Int}}`: Monomer type pairs `(i, j)` defining each bond.
- `n0i, n1i, n2i, n3i::Vector{Float64}`: Per-species scalar weighted densities.
- `nV1, nV2::SVector{3, Float64}`: Total vectorial weighted densities (directional terms).
- `σ::Vector{Float64}`: Diameters of all bead types.

# Returns
- `fᵉˣ::Float64`: Excess Helmholtz free energy density from MFMT + TPT1.

# Notes
- The MFMT term includes directional correlations via `nV1`, `nV2`.
- TPT1 chain corrections are modified by an anisotropy factor: `Ξ = 1 - |nV₂|² / n₂²`.
"""


function free_energy_TPT1mfmt(ρ_bonds, bond_types, n0i, n1i, n2i, n3i, nV1, nV2, σ)
    # Total scalar weighted densities
    n0 = sum(n0i)
    n1 = sum(n1i)
    n2 = sum(n2i)
    n3 = sum(n3i)

    # Hard-sphere contribution (includes directional correlations)
    f_hs = free_energy_mfmt(n0, n1, n2, n3, nV1, nV2)

    # Chain connectivity correction
    Xi = 1.0 - dot(nV2, nV2) / n2^2
    inv_1_minus_n3 = 1.0 / (1.0 - n3)
    term2 = n2 * Xi / 4.0 * inv_1_minus_n3^2
    term3 = n2^2 * Xi / 72.0 * inv_1_minus_n3^3

    f_ch = 0.0
    for (u, (i, j)) in enumerate(bond_types)
        σ_ij = (σ[i] + σ[j]) / 2
        τ_ij = σ[i] * σ[j] / σ_ij

        y_HS = inv_1_minus_n3 + term2 * τ_ij + term3 * τ_ij^2
        f_ch -= 0.5 * (ρ_bonds[1, u] + ρ_bonds[2, u]) * log(y_HS)
    end

    return f_hs + f_ch
end


"""
    free_energy_TPT1mfmt(rho_bonds, bond_types, n0i, n1i, n2i, n3i, σ)

Computes the excess free energy density for a bulk hard-sphere fluid using 
Modified Fundamental Measure Theory (MFMT) and TPT1 chain connectivity corrections 
without directional (vector) correlations.

# Arguments
- `rho_bonds::Vector{Float64}`: Density of each bonded pair.
- `bond_types::Vector{Tuple{Int, Int}}`: Monomer type pairs `(i, j)` defining each bond.
- `n0i, n1i, n2i, n3i::Vector{Float64}`: Per-species scalar weighted densities.
- `σ::Vector{Float64}`: Diameters of all bead types.

# Returns
- `fᵉˣ::Float64`: Excess Helmholtz free energy density from MFMT + TPT1.

# Notes
- Use this version when system is isotropic or uniform, and directional effects are negligible.
- TPT1 chain corrections depend only on `n₂`, `n₃`, and pair geometry.
"""

function free_energy_TPT1mfmt(ρ_bonds, bond_types, n0i, n1i, n2i, n3i, σ)
    # Total scalar weighted densities
    n0 = sum(n0i)
    n1 = sum(n1i)
    n2 = sum(n2i)
    n3 = sum(n3i)

    # Hard-sphere contribution from MFMT
    f_hs = free_energy_mfmt(n0, n1, n2, n3)

    # Precomputed terms for y_HS expansion
    inv_1_minus_n3 = 1.0 / (1.0 - n3)
    term2 = n2 / 4.0 * inv_1_minus_n3^2
    term3 = n2^2 / 72.0 * inv_1_minus_n3^3

    # Chain contribution
    f_ch = 0.0
    for (u, (i, j)) in enumerate(bond_types)
        σ_ij = (σ[i] + σ[j]) / 2
        τ_ij = σ[i] * σ[j] / σ_ij

        y_HS = inv_1_minus_n3 + term2 * τ_ij + term3 * τ_ij^2

        f_ch -= 0.5 * (ρ_bonds[1, u] + ρ_bonds[2, u]) * log(y_HS)
    end

    return f_hs + f_ch
end


"""
    compute_TPT1mfmt_deriv!(dΦ, lng, ρ_bonds, bond_types, n0i, n1i, n2i, n3i, nV1, nV2, σ)

Computes the MFMT + TPT1 functional derivative ∂fᵉˣ/∂nₐ for a multi-component system with chain connectivity.

# Arguments
- `dΦ::AbstractMatrix`: Output array of size (10, N) for derivatives ∂fᵉˣ/∂nₐ per species.
- `lng::AbstractVector`: Output array to store log(y_HS) values for each bond.
- `ρ_bonds::Vector{Float64}`: Bond density for each bonded pair.
- `bond_types::Vector{Tuple{Int, Int}}`: Bonded monomer type indices `(i, j)`.
- `n0i, n1i, n2i, n3i::Vector{Float64}`: Scalar weighted densities per species.
- `nV1, nV2::SVector{3, Float64}`: Total vectorial weighted densities.
- `σ::Vector{Float64}`: Diameter of each bead type.

# Notes
- The MFMT derivative is computed once using total densities and broadcast across all species.
- TPT1 corrections add contributions to `n₂`, `n₃`, and `nV₂` channels (rows 3, 4, and 8–10).
- Directionality is accounted for via the Xi prefactor: `Xi = max(1 - |nV₂|² / n₂², 0)`.

Returns the modified `dΦ` and `lng` in-place.
"""
function compute_TPT1mfmt_deriv!(
    dΦ         :: AbstractMatrix,
    lng        :: AbstractVector,
    ρ_bonds,
    bond_types,
    n0i, n1i, n2i, n3i,
    nV1, nV2,
    σ
)

    # Total weighted densities
    n0 = sum(n0i)
    n1 = sum(n1i)
    n2 = sum(n2i)
    n3 = sum(n3i)

    @. dΦ = 0.0
    @. lng = 0.0

    # Hard sphere contribution from MFMT (same across species)
    compute_mfmt_deriv!(view(dΦ, :, 1), n0, n1, n2, n3, nV1, nV2)
    @views for j in Iterators.drop(axes(dΦ, 2), 1)
        @. dΦ[:, j] = dΦ[:, 1]
    end

    # TPT1 correction terms
    Xi = max(1.0-dot(nV2, nV2)/n2^2, 0.0)
    
    inv_1_minus_n3 = 1.0 / (1.0-n3) 
    term2 = n2 * Xi / 4.0 * inv_1_minus_n3^2
    term3 = n2^2 * Xi / 72.0 * inv_1_minus_n3^3

    for (u, (i, j)) in enumerate(bond_types)
        σ_ij = (σ[i] + σ[j]) / 2
        τ_ij = σ[i] * σ[j] / σ_ij

        y_HS = inv_1_minus_n3 + term2 * τ_ij + term3 * τ_ij^2

        dlny_na2 = (τ_ij / 4.0 * inv_1_minus_n3^2 +
                        n2 * τ_ij^2 * inv_1_minus_n3^3 / 36.0) / y_HS

        dlny_na3 = (inv_1_minus_n3^2 +
                    n2 * Xi / 2.0 * inv_1_minus_n3^3 * τ_ij +
                    n2^2 * Xi / 24.0 * inv_1_minus_n3^4 * τ_ij^2) / y_HS

        for W in eachindex(n2i)
            dΦ[3, W] -= 0.5 * (ρ_bonds[1, u] + ρ_bonds[2, u]) * dlny_na2
            dΦ[4, W] -= 0.5 * (ρ_bonds[1, u] + ρ_bonds[2, u]) * dlny_na3
        end

        for v = 1:3
            dlny_naV2v = -nV2[v] * (τ_ij / (2.0 * n2) * inv_1_minus_n3^2 +
                                    τ_ij^2 / 36.0 * inv_1_minus_n3^3) / y_HS

            for W in eachindex(n2i)
                dΦ[v+7, W] -= 0.5 * (ρ_bonds[1, u] + ρ_bonds[2, u]) * dlny_naV2v
            end
        end

        lng[u] = log(y_HS)
    end

    return dΦ, lng
end


"""
    compute_TPT1mfmt_deriv!(dΦ, lng, ρ_bonds, bond_types, n0i, n1i, n2i, n3i, σ)

Computes the MFMT + TPT1 functional derivative with respect to scalar weighted densities
for a bulk (homogeneous) multi-component system, including chain connectivity corrections.

# Arguments
- `dΦ::AbstractMatrix`: Output array of size (10, N) storing ∂fᵉˣ/∂nᵅ for each species (rows 1–4).
- `lng::AbstractVector`: Output array to store log(y_HS) values for each bond.
- `rho_bonds::Vector{Float64}`: Bond density for each bonded pair.
- `bond_types::Vector{Tuple{Int, Int}}`: List of bonded monomer type pairs `(i, j)`.
- `n0i, n1i, n2i, n3i::Vector{Float64}`: Per-species scalar weighted densities.
- `σ::Vector{Float64}`: Bead diameters.

# Notes
Assumes a homogeneous system. Vectorial weighted densities are zero, and directional corrections are omitted.
"""
function compute_TPT1mfmt_deriv!(
    dΦ        :: AbstractMatrix,
    lng       :: AbstractVector,
    ρ_bonds,
    bond_types,
    n0i, n1i, n2i, n3i,
    σ
)

    # Total weighted densities
    n0 = sum(n0i)
    n1 = sum(n1i)
    n2 = sum(n2i)
    n3 = sum(n3i)

    @. dΦ = 0.0
    @. lng = 0.0

    # Hard sphere contribution from MFMT (same across species)
    compute_mfmt_deriv!(view(dΦ, :, 1), n0, n1, n2, n3)
    @views for j in Iterators.drop(axes(dΦ, 2), 1)
        @. dΦ[:, j] = dΦ[:, 1]
    end

    # TPT1 correction terms (bulk)
    Xi = 1.0  # directionality correction ≈ 1 in bulk
    
    inv_1_minus_n3 = 1.0 / (1.0-n3) 
    term2 = n2 * Xi / 4.0 * inv_1_minus_n3^2
    term3 = n2^2 * Xi / 72.0 * inv_1_minus_n3^3

    for (u, (i, j)) in enumerate(bond_types)
        σ_ij = (σ[i] + σ[j]) / 2
        τ_ij = σ[i] * σ[j] / σ_ij

        y_HS = inv_1_minus_n3 + term2 * τ_ij + term3 * τ_ij^2

        dlny_na2 = (τ_ij / 4.0 * inv_1_minus_n3^2 +
                        n2 * τ_ij^2 * inv_1_minus_n3^3 / 36.0) / y_HS

        dlny_na3 = (inv_1_minus_n3^2 +
                    n2 * Xi / 2.0 * inv_1_minus_n3^3 * τ_ij +
                    n2^2 * Xi / 24.0 * inv_1_minus_n3^4 * τ_ij^2) / y_HS

        for W in eachindex(n2i)
            dΦ[3, W] -= 0.5 * (ρ_bonds[1, u] + ρ_bonds[2, u]) * dlny_na2
            dΦ[4, W] -= 0.5 * (ρ_bonds[1, u] + ρ_bonds[2, u]) * dlny_na3
        end
        
        lng[u] = log(y_HS)
    end
end