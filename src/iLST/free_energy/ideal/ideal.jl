"""
    struct IdealFreeEnergy <: ExcessFreeEnergy

Holds reusable buffers for the ideal free energy contribution.

# Fields
- `mu::Vector{Float64}`: Reusable buffer for species-level ideal chemical potentials.
"""
struct IdealFreeEnergy <: AbstractFreeEnergy
    mu :: Vector{Float64}
end

"""
    construct_free_energy(::Type{IdealFreeEnergy}, rho::BulkDensities) -> IdealFreeEnergy

Constructs an `IdealFreeEnergy` object, initializing internal buffers
based on the number of species in the system.

# Arguments
- `rho::BulkDensities`: Bulk density structure containing `.species`.

# Returns
- `IdealFreeEnergy`: Struct with a zero-initialized chemical potential buffer.
"""
function construct_free_energy(::Type{IdealFreeEnergy}, rho::BulkDensities)
    return IdealFreeEnergy(zeros(length(rho.species)))
end


is_ideal(::IdealFreeEnergy)    = true

"""
    free_energy(molsys :: MolecularSystem, bulk :: BulkState, fe_model::IdealFreeEnergy) -> Float64

Computes the total ideal free energy based on the species densities
using the standard expression:

    Fᵢᵈ = ∑ ρᵢ (ln(ρᵢ) - 1)
"""
function free_energy(molsys :: MolecularSystem, bulk :: BulkState, fe_model::IdealFreeEnergy)
    return compute_f_ideal(bulk.rho.species)
end

"""
    chemical_potential!(molsys::MolecularSystem, bulk::BulkState, fe_model::IdealFreeEnergy)

Computes the ideal (entropic) contribution to the chemical potentials and updates the `BulkState`.

Operations:
- Computes `μ_speciesᵢ = ln(ρ_speciesᵢ)` with a floor for small densities.
- Computes segment-level `μ_segments[i][k, j] = ln(ρ_segments[i][k, j])` (state k, segment j).
- Adds the ideal species chemical potential `μ` from `fe_model` to `bulk.mu_species`.

Notes:
- Densities below 1e-100 are assigned a chemical potential of -1e20.
- Results are stored in-place in `fe_model.mu` and `bulk.mu_segments`.
"""

function chemical_potential!(molsys :: MolecularSystem, bulk :: BulkState, fe_model::IdealFreeEnergy)
    compute_mu_ideal!(fe_model.mu, bulk.rho.species)
    @. bulk.mu_species += fe_model.mu
end

"""
    compute_f_ideal(ρ::Vector{Float64}) -> Float64

Computes the ideal part of the free energy using the species densities:

    Fᵢᵈ = ∑ ρᵢ (ln(ρᵢ) - 1)

Values below 1e-100 are excluded to avoid log(0).
"""
function compute_f_ideal(ρ)
    f_ideal = 0.0
    @inbounds for (_, ρ_i) in enumerate(ρ)
        f_ideal += ρ_i > 1e-100 ? ρ_i * (log(ρ_i) - 1) : 0.0
    end
    return f_ideal
end

"""
    compute_mu_ideal!(μ_ideal::Vector, ρ::Vector)

In-place computation of the ideal (entropic) chemical potential:

    μᵢ = ln(ρᵢ)

A hard floor is applied to prevent evaluating log(0):
- If ρᵢ ≤ 1e-100, then μᵢ = -1e20.

Arguments:
- `μ_ideal`: Output vector to store chemical potentials (modified in-place).
- `ρ`: Input vector of densities.
"""

function compute_mu_ideal!(μ_ideal, ρ)
    @inbounds for (i, ρ_i) in enumerate(ρ)
        μ_ideal[i] = ρ_i > 1e-100 ? log(ρ_i) : -1e20
    end
end


"""
    compute_mu_ideal(ρ::Vector) -> Vector{Float64}

Allocating version of `compute_mu_ideal!`. Returns a new vector of:

    μᵢ = ln(ρᵢ)

A floor is applied: ρᵢ ≤ 1e-100 ⇒ μᵢ = -1e20.
"""
function compute_mu_ideal(ρ)
    μ_ideal = zeros(Float64, size(ρ))
    compute_mu_ideal!(μ_ideal, ρ)
    return μ_ideal
end
