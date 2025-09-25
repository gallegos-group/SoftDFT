function bulk_chemical_potential!(u, molsys :: MolecularSystem, bulk :: BulkState, species_model :: freely_jointed)

    gC = compute_child_propagators(u, molsys, bulk, species_model)
    gP = compute_parent_propagators(u, gC, molsys, bulk, species_model)
    evaluate_partition_function!(u, gP, gC, molsys, bulk, species_model)

    bulk.mu_species[u] -= log(bulk.Xi[u])
end