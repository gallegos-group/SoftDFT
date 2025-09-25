"""
    struct CoulombFreeEnergy <: ExcessFreeEnergy

Represents the Coulombic contribution to the Helmholtz free energy in a bulk system.

Fields:
- `dPhi::Vector{Float64}`: Per-bead excess chemical potential contributions from Coulomb interactions.
- `Psi::Vector{Float64}`: Scalar electrostatic potential (typically uniform in bulk).

This structure is used in classical DFT or liquid-state theory to represent long-range electrostatic interactions.
"""
struct CoulombFreeEnergy <: ExcessFreeEnergy
    dPhi :: Vector{Float64}
    Psi  :: Vector{Float64}
end

"""
    construct_free_energy(::Type{CoulombFreeEnergy}, rho::BulkDensities) -> CoulombFreeEnergy

Constructs an instance of `CoulombFreeEnergy` with fields initialized to zero.

- `rho`: The bulk density structure containing bead information.

Returns:
- A `CoulombFreeEnergy` object with `dPhi` sized to match the number of beads, and `Psi` initialized to a scalar value.
"""
function construct_free_energy(::Type{CoulombFreeEnergy}, rho::BulkDensities)
    return CoulombFreeEnergy(zeros(length(rho.beads)), [0.0])
end


"""
    free_energy(molsys, bulk, fe_model) -> Float64

Computes the Coulombic contribution to the excess free energy for a bulk system.

The functional form is:

    Fᵉˡ = Ψ ⋅ ∑ᵢ ρᵢ Zᵢ

# Arguments
- `molsys`: Molecular system (contains monomer valences `Zᵢ`)
- `bulk`: Bulk state (contains densities `ρᵢ` and potential `Ψ`)
- `fe_model`: Coulombic free energy structure (not used in this call, but present for consistency)

# Returns
- `Float64`: Excess coulombic free energy
"""
function free_energy(molsys :: MolecularSystem, bulk :: BulkState, fe_model :: CoulombFreeEnergy)

    Psi = bulk.Psi[1]
    rho = bulk.rho.beads
    @unpack valences = molsys.properties.monomers

    return free_energy_coul(Psi, rho, valences)
end

"""
    chemical_potential(molsys, bulk, fe_model)

Computes the Coulombic contribution to the excess chemical potential.

Each component is given by:

    μᵢᵉˡ = Ψ ⋅ Zᵢ

# Arguments
- `molsys`: Molecular system (contains valences `Zᵢ`)
- `bulk`: Bulk state (where μᵢᵉˡ is accumulated)
- `fe_model`: Stores Coulombic potential and temporary derivatives
"""
function chemical_potential!(molsys :: MolecularSystem, bulk :: BulkState, fe_model :: CoulombFreeEnergy)

    Psi = bulk.Psi[1]
    @unpack valences = molsys.properties.monomers

    compute_mu_coul!(fe_model.dPhi, Psi, valences)

    @. bulk.mu_ex += fe_model.dPhi
end


"""
    free_energy_coul(Ψ, ρ, Z) -> Float64

Evaluates the coulombic free energy:

    Fᵉˡ = Ψ ⋅ ∑ᵢ ρᵢ Zᵢ
"""
function free_energy_coul(Ψ, ρ, Z)
    return Ψ * sum(@~ @. ρ * Z)
end

"""
    compute_mu_coul!(dΦ, Ψ, Z)

Compute the Coulombic contribution to the excess chemical potential for each bead.

# Arguments
- `dΦ::Vector{Float64}`: Output vector to store chemical potentials (modified in-place).
- `Ψ::Float64`: Electrostatic potential.
- `Z::Vector{Float64}`: Valence (charge) of each bead type.

# Formula
The Coulombic chemical potential for each bead `i` is:
μᵢ^{coul} = Ψ ⋅ Zᵢ
"""

function compute_mu_coul!(dΦ, Ψ, Z)
    @. dΦ = Ψ * Z
end