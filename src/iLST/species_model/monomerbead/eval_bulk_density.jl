"""
    eval_bulk_density!(
        u::Int,
        Xi::Float64,
        struct_molsys::MolecularSystem,
        struct_bulk::BulkState
)

Populates the (pair, beads, bonds, etc.) density field for species `u`, using the current excess
chemical potential (`mu_ex`) and protonation work (`mu_prot`). No other density fields
 are modified.

Arguments:
- `u`: Species index.
- `rho1_pairs_u`: 4D array to store the pair densities [j, seg_j, i, seg_i].
- `Xi`: Partition function normalization constant.
- `struct_molsys`: Full molecular system with configuration and monomer info.
- `struct_bulk`: Contains `mu_ex` and `mu_prot`.
"""

function eval_bulk_density!(
    u,
    struct_molsys :: MolecularSystem,
    struct_bulk :: BulkState,
    species_model :: monomerbead
)

    evaluate_partition_function!(u, struct_molsys, struct_bulk, species_model)
    accumulate_densities!(u, struct_molsys, struct_bulk, species_model)
end