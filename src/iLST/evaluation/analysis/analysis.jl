# === Bulk Analysis Utilities ===
# Includes routines for computing and analyzing chemical potentials,
# excess free energies, and thermodynamic observables from bulk solutions.

include("bulk_chemical_potential.jl")
include("excess_free_energy.jl")
include("eval_pressure.jl")

"""
    summarize_bulk_thermo(molsys, bulk) -> NamedTuple

Compute a summary of key thermodynamic quantities for a bulk system,
including the chemical potentials, excess free energy, and pressure.

# Arguments
- `molsys::MolecularSystem`: Fully constructed molecular system.
- `bulk::BulkState`: Solved bulk density and field information.

# Returns
A `NamedTuple` with fields:
- `chemical_potential`: Vector of species chemical potentials
- `excess_free_energy`: Excess Helmholtz free energy (units depend on model)
- `pressure`: Bulk pressure (typically in reduced or real units)

# Example
```julia
results = summarize_bulk_thermo(system, bulk)
println(results.pressure)
"""

function summarize_bulk_thermo(molsys::MolecularSystem, bulk::BulkState)
    μ = bulk_chemical_potential(molsys, bulk)
    Fex = excess_free_energy(molsys, bulk)
    P = eval_pressure(molsys, bulk)
    return (; chemical_potential = μ, excess_free_energy = Fex, pressure = P)
end
