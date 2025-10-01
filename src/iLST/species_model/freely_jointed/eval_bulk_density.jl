"""
    eval_bulk_density(u, rho1_beads, rho1_bonds, molsys, bulk)

Evaluates the propagator-based density contributions for species `u` and
updates the bead and bond densities in-place.

This includes:
- Computing forward (`gC`) and backward (`gP`) propagators
- Evaluating the partition function `Xi`
- Updating the species chemical potential
- Accumulating density contributions

Arguments:
- `u`: Index of the species
- `rho1_beads`: Accumulator for bead densities
- `rho1_bonds`: Accumulator for bond densities
- `molsys`: Molecular system
- `bulk`: Thermodynamic state
"""
function eval_bulk_density!(
    u,
    molsys :: MolecularSystem,
    bulk :: BulkState,
    species_model :: freely_jointed
)

    gC = compute_child_propagators(u, molsys, bulk, species_model)
    gP = compute_parent_propagators(u, gC, molsys, bulk, species_model)
    evaluate_partition_function!(u, gP, gC, molsys, bulk, species_model)
    accumulate_densities!(u, gP, gC, molsys, bulk, species_model)
end