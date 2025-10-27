include("excess_free_energy.jl")
"""
    calculate_bulk_pressure(molsys, bulk::BulkState)

Computes the excess free energy and updates the bulk pressure field
in-place within `bulk`.

This function does not return anything.
"""
function eval_bulk_pressure(molsys :: MolecularSystem, bulk :: BulkState)

    rho_species = bulk.rho.species
    rho_beads   = bulk.rho.beads
    rho_bonds   = bulk.rho.bonds

    mu_ex = bulk.mu_ex
    lng   = bulk.lng

    f_ex = 0.0 # excess_free_energy(molsys, bulk)

    bulk.Pressure[1] = sum(rho_species) - f_ex +
                              sum(@. @~ rho_beads * mu_ex) -
                              sum(@views @. @~ rho_bonds[1, :] * lng) / 2.0 -
                              sum(@views @. @~ rho_bonds[2, :] * lng) / 2.0
end