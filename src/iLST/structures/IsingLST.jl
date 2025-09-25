"""
    IsingLST

Encapsulates all data required to solve a liquid-state theory (LST) problem
for weak polyelectrolytes using Ising Density Functional Theory.

# Fields
- `molsys::MolecularSystem`: Full specification of species, reactions, and configurations.
- `bulk::BulkState`: Bulk thermodynamic state, densities, and potentials.
- `numerics::Dict{String, Float64}`: Numerical solver parameters (e.g., damping, tolerance).
"""

include("BulkState.jl")

struct IsingLST #<: Procedure_iLST
    molsys     :: MolecularSystem               # Molecular system specification
    bulk       :: BulkState                     # Bulk thermodynamic problem and state
    numerics   :: Dict{String, Float64}         # Input numerical parameters
end
