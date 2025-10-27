
"""
    eval_bulk_density(molsys, bulk) -> (rho1_beads, rho1_bonds)

Computes updated bead and bond densities for the bulk system by looping
over all species and dispatching to their corresponding model-specific
density evaluation routines.

Raises:
- `ErrorException` if any NaNs are detected in the resulting densities.
"""
function eval_bulk_density!(
    molsys::MolecularSystem,
    bulk::BulkState
)
 
    @. bulk.rho.beads = 0.0
    @. bulk.rho.bonds = 0.0
    
    for (u, config) in enumerate(molsys.configurations)
        species_model = config.species_model

        eval_bulk_density!(u, molsys, bulk, species_model)
    end
end

include("../../../species_model/analyticalEval.jl")

"""
    eval_bulk_density(::Int, ...) for unimplemented species model

Fallback method that throws an error if no implementation exists
for the provided `species_model`.
"""
function eval_bulk_density!(
    u,
    molsys::MolecularSystem,
    bulk::BulkState,
    species_model
)
    error("No `eval_bulk_density` implementation for species model of type `$(typeof(species_model))`.")
end
