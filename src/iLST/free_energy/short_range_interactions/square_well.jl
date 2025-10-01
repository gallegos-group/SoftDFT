"""
    struct SquareWellFreeEnergy <: ExcessFreeEnergy

Holds the excess chemical potential contributions from square-well interactions.

# Fields
- `dPhi::Vector{Float64}`: Functional derivative ∂f/∂ρᵢ for each bead type.
"""
struct SquareWellFreeEnergy <: ExcessFreeEnergy
    dPhi :: Vector{Float64}
end

"""
    construct_free_energy(::Type{SquareWellFreeEnergy}, rho::BulkDensities) -> SquareWellFreeEnergy

Initializes the `SquareWellFreeEnergy` structure with zeroed chemical potential derivatives.

# Arguments
- `rho::BulkDensities`: Contains bead density arrays used to size `dPhi`.

# Returns
- `SquareWellFreeEnergy`: Instance with `dPhi` initialized to zeros.
"""
function construct_free_energy(::Type{SquareWellFreeEnergy}, rho::BulkDensities)
    return SquareWellFreeEnergy(zeros(Float64, length(rho.beads)))
end



"""
    free_energy(molsys::MolecularSystem, bulk::BulkState, fe_model::SquareWellFreeEnergy) -> Float64

Computes the square-well excess free energy using mean-field theory.

# Arguments
- `molsys`: Molecular system with bead diameters and pair interaction parameters.
- `bulk`: Bulk state containing bead densities.
- `fe_model`: Square-well interaction structure (not used here directly).

# Returns
- Scalar free energy contribution from square-well attractions.
"""

function free_energy(molsys::MolecularSystem, bulk::BulkState, fe_model :: SquareWellFreeEnergy)
    rho = bulk.rho.beads

    @unpack properties = molsys
    @unpack diameters = properties.monomers
    @unpack species, pairs, lambdas, energys = properties.pairs[:square_well]

    return free_energy_sw(rho, diameters, species, pairs, lambdas, energys)
end

"""
    free_energy_sw(ρ, σ, sw_species, sw_pairs, λ_ij, ε_ij) -> Float64

Computes the mean-field square-well excess free energy.

# Arguments
- `ρ`: Vector of bead densities.
- `σ`: Vector of bead diameters.
- `sw_species`: Indices of species participating in square-well interactions.
- `sw_pairs`: Tuple pairs indicating interacting species.
- `λ_ij`: Well width multipliers.
- `ε_ij`: Well depth (energy) values.

# Returns
- Scalar excess free energy from square-well attractions.
"""

function free_energy_sw(ρ, σ, sw_species, sw_pairs, λ_ij, ε_ij)
    f_sw = 0.0
    for i in sw_species
        for j in sw_species
            index = findfirst(x -> x == (min(i, j), max(i, j)), sw_pairs)
            isnothing(index) && continue

            σ_ij = (σ[i] + σ[j]) / 2.0
            V_ij = (4.0*π/3.0) * σ_ij^3 * (λ_ij[index]^3 - 1.0)

            f_sw -= 0.5 * ρ[i] * ρ[j] * ε_ij[index] * V_ij
        end
    end

    return f_sw
end


"""
    chemical_potential!(molsys, bulk, fe_model)

Computes the functional derivative of the square-well free energy and updates `mu_ex`.

# Arguments
- `molsys`: Molecular system with monomer properties.
- `bulk`: Bulk state (provides densities and `mu_ex` to be updated).
- `fe_model`: Square-well structure that stores the computed derivative.

- Updates `bulk.mu_ex` in-place by adding the square-well contribution.
"""

function chemical_potential!(molsys::MolecularSystem, bulk::BulkState, fe_model :: SquareWellFreeEnergy)

    rho = bulk.rho.beads

    @unpack properties = molsys
    @unpack diameters = properties.monomers
    @unpack species, pairs, lambdas, energys = properties.pairs[:square_well]
    @unpack dPhi = fe_model

    compute_sw_deriv!(dPhi, rho, diameters, species, pairs, lambdas, energys)
    
    @. bulk.mu_ex += dPhi
end

"""
    compute_sw_deriv!(dΦ, ρ, σ, sw_species, sw_pairs, λ_ij, ε_ij)

Computes the functional derivative ∂f/∂ρᵢ of the square-well free energy and stores in `dΦ`.

# Arguments
- `dΦ`: Output vector of chemical potential contributions (updated in-place).
- `ρ`: Vector of bead densities.
- `σ`: Vector of bead diameters.
- `sw_species`: Indices of species in square-well interactions.
- `sw_pairs`: Vector of tuple pairs (i, j) of species that interact.
- `λ_ij`: Width multipliers for each interaction pair.
- `ε_ij`: Well depth (energy) for each pair.

# Returns
- The updated `dΦ` vector.
"""

function compute_sw_deriv!(dΦ, ρ, σ, sw_species, sw_pairs, λ_ij, ε_ij)
    @. dΦ = 0.0

    for i in sw_species
        for j in sw_species
            index = findfirst(x -> x == (min(i, j), max(i, j)), sw_pairs)
            isnothing(index) && continue

            σ_ij = (σ[i] + σ[j]) / 2.0
            V_ij = (4.0*π / 3.0) * σ_ij^3 * (λ_ij[index]^3 - 1.0)

            dΦ[i] -= ρ[j] * ε_ij[index] * V_ij
        end
    end

    return dΦ
end
