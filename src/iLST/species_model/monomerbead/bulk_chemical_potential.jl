function bulk_chemical_potential!(u, struct_molsys :: MolecularSystem, struct_bulk :: BulkState, species_model :: monomerbead)

    evaluate_partition_function!(u, struct_molsys, struct_bulk, species_model)

    struct_bulk.mu_species[u] -= log(struct_bulk.Xi[u])
end