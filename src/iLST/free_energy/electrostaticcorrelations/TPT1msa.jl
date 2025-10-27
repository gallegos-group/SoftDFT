"""
    struct TPT1MSAFreeEnergy <: ExcessFreeEnergy

Represents the TPT1-enhanced Mean Spherical Approximation (MSA) contribution to the excess free energy.

# Fields
- `dPhi::Vector{Float64}`: Per-bead excess chemical potential.
- `lng::Vector{Float64}`: Logarithmic bonding contributions (TPT1 correction).
- `wt_rho::Vector{Float64}`: Working array for weighted bead densities.
- `wt_bonds::Matrix{Float64}`: Working array for weighted bond densities.
"""
struct TPT1MSAFreeEnergy <: ExcessFreeEnergy
    dPhi     :: Vector{Float64}
    lng      :: Vector{Float64}
    wt_rho   :: Vector{Float64}
    wt_bonds :: Matrix{Float64}
end

"""
    construct_free_energy(::Type{TPT1MSAFreeEnergy}, rho::BulkDensities) -> TPT1MSAFreeEnergy

Initializes a `TPT1MSAFreeEnergy` instance with zeroed fields based on the bead and bond counts.

# Arguments
- `rho::BulkDensities`: Bulk density object containing bead and bond information.

# Returns
- A `TPT1MSAFreeEnergy` object with all fields initialized to zero.
"""
function construct_free_energy(::Type{TPT1MSAFreeEnergy}, rho::BulkDensities)
    nbeads = length(rho.beads)
    nbonds = size(rho.bonds, 2)

    return TPT1MSAFreeEnergy(zeros(nbeads), zeros(nbonds), zeros(nbeads), zeros(2, nbonds))
end


"""
    free_energy(molsys::MolecularSystem, bulk::BulkState, fe_model::TPT1MSAFreeEnergy) -> Float64

Computes the excess Helmholtz free energy contribution from the TPT1-corrected MSA (Mean Spherical Approximation) model.

# Arguments
- `molsys::MolecularSystem`: Molecular system structure containing monomer and bond properties.
- `bulk::BulkState`: Bulk state containing bead and bond densities.
- `fe_model::TPT1MSAFreeEnergy`: Structure holding intermediate values for the TPT1-MSA model.

# Returns
- `Float64`: Total excess free energy contribution from the TPT1-MSA model.
"""

function free_energy(molsys::MolecularSystem, bulk::BulkState, fe_model :: TPT1MSAFreeEnergy)
    # Densities
    rho = bulk.rho.beads
    rho_bonds = bulk.rho.bonds

    # Properties
    @unpack properties = molsys
    @unpack bond_types = properties
    @unpack diameters, valences = properties.monomers
    @unpack bjerrum_length = properties.system

    return free_energy_TPT1msa(rho_bonds, bond_types, rho, diameters, valences, bjerrum_length)
end

"""
    chemical_potential!(molsys::MolecularSystem, bulk::BulkState, fe_model::TPT1MSAFreeEnergy)

Computes the excess chemical potentials and logarithmic bonding terms from the TPT1-corrected MSA model, and stores them in the bulk state.

# Arguments
- `molsys::MolecularSystem`: Molecular system structure containing monomer and bond properties.
- `bulk::BulkState`: Bulk state containing densities and storage arrays for chemical potentials.
- `fe_model::TPT1MSAFreeEnergy`: Structure holding working arrays and output buffers (`dPhi`, `lng`).

# Side Effects
- Updates `bulk.mu_ex` with excess chemical potentials.
- Updates `bulk.lng` with logarithmic bonding corrections.
"""

function chemical_potential!(molsys::MolecularSystem, bulk::BulkState, fe_model :: TPT1MSAFreeEnergy)
    # Densities
    rho = bulk.rho.beads
    rho_bonds = bulk.rho.bonds

    # Properties
    @unpack properties = molsys
    @unpack bond_types = properties
    @unpack diameters, valences = properties.monomers
    @unpack bjerrum_length = properties.system
    @unpack dPhi, lng = fe_model

    compute_TPT1msa_deriv!(dPhi, lng, rho_bonds, bond_types, rho, diameters, valences, bjerrum_length)

    @. bulk.mu_ex += dPhi
    @. bulk.lng += lng
end

"""
    free_energy_TPT1msa(ρ_bonds, bond_types, ρ, σ, Z, ℓB) -> Float64

Computes the total excess free energy for an ionic chain fluid using 
TPT1 corrections to the Mean Spherical Approximation (MSA).

# Arguments
- `ρ_bonds`: Bond densities for each bonded pair.
- `bond_types`: Vector of index pairs `(i, j)` representing bonded bead types.
- `ρ`: Bead densities for each species.
- `σ`: Bead diameters.
- `Z`: Valences (charges) of each bead type.
- `ℓB`: Bjerrum length.

# Returns
- Total excess free energy density (in units of kBT per unit volume).

# Notes
The total excess free energy is given by:

    Fᵉˣ = Fᵉˡ + Fᶜʰ

where:

    Fᵉˡ = Γ³ / (3π) - ℓB ⋅ ∑ᵢ ρᵢ Zᵢ (Zᵢ Γ + η σᵢ) / (1 + Γ σᵢ)

and

    Fᶜʰ = -½ ⋅ ∑ᵤ ρᵦᵤ ⋅ ln gᵤ

with

    ln gᵤ = -ℓB / Dᵢⱼ ⋅ Zₑᶠᶠ(i) ⋅ Zₑᶠᶠ(j)

    Zₑᶠᶠ(i) = (Zᵢ - η σᵢ²) / (1 + Γ σᵢ)

    σᵢⱼ = (σᵢ + σⱼ) / 2
"""

function free_energy_TPT1msa(ρ_bonds, bond_types, ρ, σ, Z, ℓB)
    f_el = free_energy_msa(ρ, σ, Z, ℓB)
    Γ, η, _ = Gamma_msa(ρ, σ, Z, ℓB)

    f_ch = 0.0
    for (u, (i, j)) in enumerate(bond_types)
        σ_ij = (σ[i] + σ[j]) / 2.0
        
        X_i  = (Z[i] - η*σ[i]^2) / (1.0 + Γ*σ[i])
        X_j  = (Z[j] - η*σ[j]^2) / (1.0 + Γ*σ[j])

        lng = -ℓB / σ_ij * X_i * X_j

        f_ch -= 0.5 * (ρ_bonds[1, u] + ρ_bonds[2, u]) * lng
    end

    return f_el + f_ch
end

"""
    compute_TPT1msa_deriv!(dΦ, lng, ρ_bonds, bond_types, ρ, σ, Z, ℓB)

Computes the functional derivative of the TPT1-MSA excess free energy with respect to bead densities.

# Arguments
- `dΦ::Vector{Float64}`: Output vector for chemical potential contributions ∂Fᵉˣ/∂ρᵢ (modified in-place).
- `lng::Vector{Float64}`: Output vector for `log(gᵢⱼ)` values from electrostatic correction (modified in-place).
- `ρ_bonds::Vector{Float64}`: Bond densities for each bonded pair.
- `bond_types::Vector{Tuple{Int, Int}}`: List of bonded bead type pairs `(i, j)`.
- `ρ::Vector{Float64}`: Bead number densities.
- `σ::Vector{Float64}`: Bead diameters.
- `Z::Vector{Float64}`: Bead valences (charges).
- `ℓB::Float64`: Bjerrum length.

# Description
This function computes:

1. The MSA free energy derivative via `compute_msa_deriv!`.
2. The chain connectivity correction via the electrostatic TPT1 approximation.

The electrostatic contribution is:

    ln gᵢⱼ = -ℓB / σᵢⱼ ⋅ X(i) ⋅ X(j)

with:

    X(i) = (Zᵢ - η σᵢ²) / (1 + Γ σᵢ)

where `Γ` and `η` are the MSA screening parameters.

Their implicit ρ-dependence is handled using chain rule derivatives:

- ∂X/∂ρ accounts for ∂Γ/∂ρ and ∂η/∂ρ,
- ∂Γ/∂ρ is computed from its self-consistent closure relation.

The correction is accumulated in `dΦ` as:

    dΦᵢ += ∑ⱼ ½ ⋅ ρ_bond ⋅ ∂ln gⱼ/∂ρᵢ

Returns `dΦ` and `lng`, both modified in-place.
"""

function compute_TPT1msa_deriv!(dΦ, lng, ρ_bonds, bond_types, ρ, σ, Z, ℓB)

    @. lng = 0.0
    @. dΦ = 0.0

    compute_msa_deriv!(dΦ, ρ, σ, Z, ℓB)

    Γ, η, H = Gamma_msa(ρ, σ, Z, ℓB)

    if isapprox(Γ, 0.0)
        for (u, (i, j)) in enumerate(bond_types)
            σ_ij = (σ[i] + σ[j]) / 2.0

            X_i = Z[i]
            X_j = Z[j]

            lng[u] = -ℓB / σ_ij * X_i * X_j
        end

        return
    else
        for (u, (i, j)) in enumerate(bond_types)
            σ_ij = (σ[i] + σ[j]) / 2.0

            X_i = (Z[i] - η*σ[i]^2) / (1.0 + Γ*σ[i])
            X_j = (Z[j] - η*σ[j]^2) / (1.0 + Γ*σ[j])

            lng[u] = -ℓB / σ_ij * X_i * X_j
        end

        ∂f_Γ = -2.0*π*ℓB * sum(@~ @. ρ*σ * 
                    (Z - η*σ*σ)*(Z - η*σ*σ)/(1.0 + Γ*σ)/(1.0 + Γ*σ)/(1.0 + Γ*σ))

        ∂f_η = -2.0*π*ℓB * sum(@~ @. ρ*σ*σ * (Z - η*σ*σ)/(1.0 + Γ*σ)/(1.0 + Γ*σ))

        ∂η_Γ = -1.0/H * sum(@~ @. ρ*Z*σ*σ / (1.0 + Γ*σ)/(1.0 + Γ*σ))

        ∂η_H = -η / H

        ∂H_Γ = -sum(@~ @. ρ*σ*σ*σ*σ / (1.0 + Γ*σ) / (1.0 + Γ*σ))

        dη_Γ = ∂η_Γ + ∂η_H*∂H_Γ

        denom = 1.0 / (2.0*Γ - ∂f_η*dη_Γ - ∂f_Γ)
        for w in eachindex(ρ)
            inv_1_plus_Γσ = 1.0 / (1.0 + Γ*σ[w])

            X_w = (Z[w] - η*σ[w]^2) * inv_1_plus_Γσ

            ∂f_ρ = π*ℓB * X_w^2
            ∂η_ρ = σ[w] / H * Z[w] * inv_1_plus_Γσ
            ∂H_ρ = σ[w]^3 * (inv_1_plus_Γσ - 1.0/3.0)

            dΓ_ρ = (∂f_ρ + ∂f_η * (∂η_ρ + ∂f_η*∂η_H*∂H_ρ)) * denom
            dη_ρ = dη_Γ * dΓ_ρ

            for (u, (i, j)) in enumerate(bond_types)
                σ_ij = (σ[i] + σ[j]) / 2.0
                
                # J terms
                inv_1_plus_Γσi = 1.0/(1.0+Γ*σ[i])
                X_i = (Z[i]-η*σ[i]^2)*inv_1_plus_Γσi

                ∂Xi_η = σ[i]^2 * inv_1_plus_Γσi
                ∂Xi_Γ = σ[i] * X_i * inv_1_plus_Γσi

                dXi_ρ = -(∂Xi_η * dη_ρ + ∂Xi_Γ * dΓ_ρ)

                # J1 terms
                inv_1_plus_Γσj = 1.0 / (1.0 + Γ*σ[j]) 
                X_j = (Z[j] - η*σ[j]^2)*inv_1_plus_Γσj

                ∂Xj_η = σ[j]^2 * inv_1_plus_Γσj
                ∂Xj_Γ = σ[j] * X_j * inv_1_plus_Γσj

                dXj_ρ = -(∂Xj_η * dη_ρ + ∂Xj_Γ * dΓ_ρ)

                dΦ[w] += (ρ_bonds[1, u] + ρ_bonds[2, u]) / 2.0 * 
                            ℓB / σ_ij * (dXi_ρ*X_j + X_i*dXj_ρ)
            end
        end
    end
end