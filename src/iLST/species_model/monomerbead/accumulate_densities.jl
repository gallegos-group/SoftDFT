"""
    accumulate_densities!(
        u, Xi, struct_molsys, struct_bulk
    )

Accumulates self-pair segment–state densities for species `u` modeled as a
non-interacting monomer with discrete internal states.

This routine computes:
- Diagonal pair density entries ρ₁(i,i) for each internal state of the monomer,
  normalized by the total partition function `Xi`.
- These are stored in the 4D array `rho1_pairs_u`, where only entries of the form
  `rho1_pairs_u[idx, seg, idx, seg]` are modified.

Arguments:
- `u`: Index of the species being evaluated.
- `rho1_pairs_u`: 4D array to accumulate state–segment pair densities (modified in-place).
- `Xi`: Partition function normalization constant.
- `struct_molsys`: Molecular system object containing species configurations and parameters.
- `struct_bulk`: Bulk thermodynamic state, including excess chemical potentials and densities.

Notes:
- This version assumes an unconnected (non-polymeric) species, such as a single monomer with ionizable states.
- Off-diagonal terms in `rho1_pairs_u` are left untouched (set to zero).
- Bond and segment connectivity are ignored.
"""



function accumulate_densities!(
    u::Int,
    struct_molsys::MolecularSystem,
    struct_bulk::BulkState,
    species_model::monomerbead
)
    config = struct_molsys.configurations[u]
    @unpack sequence, state_family = config
    @unpack mu_ex, rho = struct_bulk
    @unpack delta_muH = struct_molsys.properties.monomers

    Xi = struct_bulk.Xi[u]

    rho_pairs_u = rho.pairs[u]
    rho_segments_u = rho.segments[u]
    rho_beads = rho.beads

    n_segments = length(sequence)

    @. rho_pairs_u = 0.0  # zero existing values
    @. rho_segments_u = 0.0  # zero existing values

    for seg_i in 1:n_segments
        for (idx_i, state_i) in enumerate(state_family[seg_i])
            weight_i = exp(-mu_ex[state_i] - delta_muH[state_i])

            val = (rho.species[u] / Xi) * weight_i

            rho_pairs_u[idx_i, seg_i, idx_i, seg_i] = val
            
            rho_beads[state_i] += val
            rho_segments_u[idx_i, seg_i] += val
        end
    end
end