"""
    struct MSAFreeEnergy <: ExcessFreeEnergy

Represents the Mean Spherical Approximation (MSA) contribution to the excess free energy.

# Fields
- `dPhi::Vector{Float64}`: Per-bead excess chemical potential from MSA.
- `wt_rho::Vector{Float64}`: Working array for storing weighted bead densities (can be reused across evaluations).
"""
struct MSAFreeEnergy <: ExcessFreeEnergy
    dPhi   :: Vector{Float64}
    wt_rho :: Vector{Float64}
end

"""
    construct_free_energy(::Type{MSAFreeEnergy}, rho::BulkDensities) -> MSAFreeEnergy

Initializes an `MSAFreeEnergy` instance with zeroed fields based on the number of bead types.

# Arguments
- `rho::BulkDensities`: Bulk density structure containing bead density information.

# Returns
- A `MSAFreeEnergy` object with zero-initialized fields.
"""
function construct_free_energy(::Type{MSAFreeEnergy}, rho::BulkDensities)
    nbeads = length(rho.beads)
    return MSAFreeEnergy(zeros(nbeads), zeros(nbeads))
end


"""
    free_energy(molsys, bulk, fe_model) -> Float64

Computes the MSA contribution to the excess free energy for a bulk system.

The form used is:
    Fᶜᵒᵘˡ = analytic expression from MSA theory

# Arguments
- `molsys`: Molecular system with `diameters`, `valences`, and `bjerrum_length`.
- `bulk`: Bulk state with bead densities.
- `fe_model`: MSA structure (not directly used here but passed for consistency).

# Returns
- Total excess free energy from MSA.
"""
function free_energy(molsys::MolecularSystem, bulk::BulkState, fe_model::MSAFreeEnergy)
    # Bead density
    rho = bulk.rho.beads

    # Properties
    @unpack properties = molsys
    @unpack diameters, valences = properties.monomers
    @unpack bjerrum_length = properties.system

    return free_energy_msa(rho, diameters, valences, bjerrum_length)
end

"""
    chemical_potential!(molsys, bulk, fe_model)

Computes the MSA contribution to the excess chemical potential for each bead
and accumulates it into `mu_ex`.

# Arguments
- `molsys`: Molecular system with valences and diameters.
- `bulk`: Contains densities and the `mu_ex` field to be updated.
- `fe_model`: Structure storing the output `dPhi`.

# Notes
Modifies `bulk.mu_ex` in-place.
"""
function chemical_potential!(molsys::MolecularSystem, bulk::BulkState, fe_model::MSAFreeEnergy)
    # Bead density
    rho = bulk.rho.beads

    # Properties
    @unpack properties = molsys
    @unpack diameters, valences = properties.monomers
    @unpack bjerrum_length = properties.system

    @unpack dPhi = fe_model

    compute_msa_deriv!(dPhi, rho, diameters, valences, bjerrum_length)

    @. bulk.mu_ex += dPhi
end

"""
    free_energy_msa(ρ, σ, Z, ℓB) -> Float64

Computes the excess electrostatic free energy density of an ionic fluid using 
the Mean Spherical Approximation (MSA).

# Arguments
- `ρ`: Vector of species densities.
- `σ`: Vector of hard-sphere diameters for each species.
- `Z`: Vector of valences (charges) for each species.
- `ℓB`: Bjerrum length (in units consistent with `σ`).

# Returns
- Excess electrostatic free energy density (in units of kBT per unit volume).

# Notes
The expression used is:

    F_coul = Γ³ / (3π) - ℓB × ∑ᵢ ρᵢ Zᵢ (Zᵢ Γ + η σᵢ) / (1 + Γ σᵢ)

where Γ and η are MSA screening parameters determined self-consistently
from the system state via `Gamma_msa`.
"""


function free_energy_msa(ρ, σ, Z, ℓB)
    Γ, η, _ = Gamma_msa(ρ, σ, Z, ℓB)
    return Γ^3 / (3.0*π) - ℓB * sum(@~ @. ρ*Z * (Z*Γ + η*σ) / (1.0 + Γ*σ))
end

"""
    compute_msa_deriv!(dΦ, ρ, σ, Z, ℓB)

Computes the electrostatic excess chemical potential for each species using the 
Mean Spherical Approximation (MSA), including charge renormalization effects.

# Arguments
- `dΦ::Vector{Float64}`: Output vector of excess chemical potentials (modified in-place).
- `ρ::Vector{Float64}`: Bulk densities of each species.
- `σ::Vector{Float64}`: Hard-sphere diameters.
- `Z::Vector{Float64}`: Valences (charges) of the species.
- `ℓB::Float64`: Bjerrum length.

# Formula
The chemical potential for species i is computed as:

    μᵢ = -ℓB * [ Γ Zᵢ² / (1 + Γ σᵢ)
               + η σᵢ ( (2 Zᵢ - η σᵢ²) / (1 + Γ σᵢ) + η σᵢ² / 3 )
               + Zᵢ u* ]

where the background energy term u* is:

    u* = -π / 6 * ∑ⱼ ρⱼ σⱼ² ( (Zⱼ - η σⱼ²) / (1 + Γ σⱼ) + Zⱼ )

Here, Γ and η are computed from `Gamma_msa`.
"""

function compute_msa_deriv!(dΦ, ρ, σ, Z, ℓB)
    Γ, η, _ = Gamma_msa(ρ, σ, Z, ℓB)

    u_star = -π/6 * sum(@~ @. ρ*σ^2 * ((Z - η*σ*σ) / (1.0 + Γ*σ) + Z))

    @. dΦ = -ℓB * (Γ*Z^2 / (1.0 + Γ*σ)
                    + η*σ * ((2.0*Z-η*σ^2) / (1.0 + Γ*σ) +η*σ^2 / 3.0) 
                    + Z*u_star)
end

"""
    Gamma_msa(ρ, σ, Z, ℓB) -> (Γ, η, H)

Computes the MSA screening parameter `Γ`, the correlation coefficient `η`, and the normalization factor `H`
for a charged hard-sphere mixture using the Mean Spherical Approximation (MSA) with charge renormalization.

# Arguments
- `ρ::Vector{Float64}`: Number density of each species.
- `σ::Vector{Float64}`: Diameter of each species.
- `Z::Vector{Float64}`: Valence (charge) of each species.
- `ℓB::Float64`: Bjerrum length.

# Returns
- `Γ::Float64`: MSA screening parameter.
- `η::Float64`: Renormalization factor for effective charge.
- `H::Float64`: Auxiliary normalization factor used in `η` computation.

# Notes
- The algorithm first computes the Debye screening parameter `κ`.
- If `κ` is below tolerance, an analytical limit is used.
- Otherwise, the screening parameter `Γ` is iteratively solved using fixed-point updates.
- `σ_eff` is the average diameter weighted by species densities.
"""

function Gamma_msa(ρ, σ, Z, ℓB)

    tol = 1e-10
    Γ1 = 0.0
    κ = sqrt(4.0*π*ℓB * sum(@~ @. ρ * Z^2))

    if κ < tol
        Γ   = κ / 2.0
        η   = 0.0
        H   = 2.0 / π
    else
        σ_eff = sum(@~ @. ρ*σ)/sum(ρ)

        x = σ_eff * sqrt(4.0*π*ℓB * sum(@~ @. ρ*Z*Z))
        Γ = (sqrt(1.0 + 2.0*x) - 1.0) / (2.0*σ_eff)

        term2 = 2.0/pi*(1.0 - π/6.0 * sum(@~ @. ρ*σ*σ*σ))
    
        max_iter = 1000
        iter = 0

        while true
            H = sum(@~ @. ρ * σ^3 / (1.0 + Γ * σ)) + term2

            η = sum(@~ @. ρ * Z * σ / (1.0 + Γ * σ)) / H

            Γ1 = sqrt(π * ℓB * sum(@~ @. ρ * (Z - η * σ^2)^2 / (1.0 + Γ * σ)^2))

            # Update Γ with damping
            Γ_new = 0.9 * Γ + 0.1 * Γ1

            iter += 1
            if abs(Γ1 - Γ) < tol || iter ≥ max_iter
                Γ = Γ_new
                break
            end

            Γ = Γ_new
        end

        if iter == max_iter
            if isnan(abs(Γ1 - Γ))
                error("Γ is NaN")
            end
            @warn "Γ iteration reached max_iter = $max_iter without convergence (Δ = $(abs(Γ1 - Γ)))"
        end
    end

    return Γ, η, H
end