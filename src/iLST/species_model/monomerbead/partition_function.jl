"""
    evaluate_partition_function!(
        u::Int,
        struct_molsys::MolecularSystem,
        bulk::BulkState,
        species_model::monomerbead,
    ) -> Xi

Computes the partition function `Xi` for a species modeled as a single monomer bead.
This is used to normalize the density contribution for simple species without internal structure.

Arguments:
- `u`: Index of the species/configuration.
- `struct_molsys`: Molecular system description, including sequence and protonation energies.
- `bulk`: Bulk thermodynamic state containing excess chemical potentials.
- `species_model`: Dispatch type `monomerbead` (used for method selection).
"""

function evaluate_partition_function!(
    u::Int,
    struct_molsys::MolecularSystem,
    bulk::BulkState,
    species_model :: monomerbead
)
    @unpack configurations, properties = struct_molsys
    @unpack mu_ex = bulk
    @unpack delta_muH = properties.monomers

    config = configurations[u]
    @unpack sequence, state_family = config

    bulk.Xi[u] = 0.0
    root = 1  # Assumes single-monomer species starts at segment 1

    for state in state_family[root]
        bulk.Xi[u] += exp(-mu_ex[state] - delta_muH[state])
    end
end
