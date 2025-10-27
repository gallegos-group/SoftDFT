"""
    accumulate_densities!(
        u, gP, gC,
        molsys, bulk
    )

Compute the state-resolved pair densities `rho1_pairs_u` for species `u` in a bulk system.

This routine evaluates:
1. The diagonal and bonded neighbor contributions using the propagators `gP`, `gC`.
2. Off-diagonal non-neighbor terms via a mean-field estimate (for closure).

Arguments:
- `u`: Index of the species being evaluated.
- `rho1_pairs_u`: 4D array to accumulate the segment–segment pair densities.
- `gP`: Parent propagator array (3D).
- `gC`: Child propagator array (4D).
- `molsys`: Molecular system with configurations and parameters.
- `bulk`: BulkState with chemical potentials and densities.

Notes:
- Only the first child is used when computing pair propagator coupling.
- Symmetry `rho1_pairs_u[i,j,k,l] = rho1_pairs_u[k,l,i,j]` is enforced.
"""

function accumulate_densities!(
    u::Int,
    gP::Array{Float64, 3},
    gC::Array{Float64, 4},
    molsys::MolecularSystem,
    bulk::BulkState,
    species_model :: freely_jointed
)
    @unpack configurations, properties = molsys
    @unpack sequence, state_family, topology = configurations[u]
    @unpack children, parents = topology
    @unpack mu_ex, rho = bulk
    @unpack delta_muH = properties.monomers
    @unpack bond_types = properties

    Xi = bulk.Xi[u]

    n_segments = length(sequence)

    rho_pairs_u = rho.pairs[u]
    rho_segments_u = rho.segments[u]
    rho_beads = rho.beads
    rho_bonds = rho.bonds

    @. rho_pairs_u = 0.0
    @. rho_segments_u = 0.0

    for seg_i in 1:n_segments
        for (idx_i, state_i) in enumerate(state_family[seg_i])
            λ_i = mu_ex[state_i] + delta_muH[state_i]
            prefactor = rho.species[u] * exp(-λ_i) * sum(@views gP[:, idx_i, seg_i]) / Xi

            if !isempty(children[seg_i])
                child = children[seg_i][1]
               
                child_contr = 0.0
                for (idx_j, _) in enumerate(state_family[child])
                    temp = gC[idx_j, 1, idx_i, seg_i]

                    for (idx_oc, other_child) in enumerate(children[seg_i])
                        if other_child != child
                            temp *= sum(@views gC[:, idx_oc, idx_i, seg_i])
                        end
                    end
                    child_contr += temp
                end
            else
                child_contr = 1.0
            end

            val = prefactor * child_contr

            # Treat the diagonal term
            rho_pairs_u[idx_i, seg_i, idx_i, seg_i] = val
            rho_segments_u[idx_i, seg_i] = val
            rho_beads[state_i] += val

            for seg_j in 1:n_segments
                for (idx_j, state_j) in enumerate(state_family[seg_j])
                    
                    idx_c = findfirst(==(seg_j), children[seg_i])
                    is_child = idx_c !== nothing # true if found

                    if is_child
                        Q_c = sum(@views gC[:, idx_c, idx_i, seg_i])
                        f = val * gC[idx_j, idx_c, idx_i, seg_i] / Q_c
                        
                        rho_pairs_u[idx_j, seg_j, idx_i, seg_i] = f
                        rho_pairs_u[idx_i, seg_i, idx_j, seg_j] = f

                        bond = (min(state_i, state_j), max(state_i, state_j))
                        idx = findfirst(==(bond), bond_types)
                        if idx !== nothing
                            @views @. rho_bonds[:, idx] += f
                        end
                    end
                end
            end
        end
    end

    # === Mean-field estimate for off-diagonal, non-neighbor terms ===
    # ===  Can be improved by evaluating extended propagator fxns  ===
    for seg_i in 1:n_segments
        for (idx_i, _) in enumerate(state_family[seg_i])
            # Treat the non-(neighbor and diagonal) term
            rho_si = rho_pairs_u[idx_i, seg_i, idx_i, seg_i]

            for seg_j in 1:n_segments
                if seg_j == seg_i || seg_j in children[seg_i] || seg_j in parents[seg_i]
                    continue
                end

                for (idx_j, _) in enumerate(state_family[seg_j])
                    rho_j = sum(@views rho_pairs_u[:, seg_j, idx_j, seg_j])

                    if rho_j < 1e-100
                        continue
                    else
                        alpha_sj = rho_pairs_u[idx_j, seg_j, idx_j, seg_j] / rho_j

                        rho_pairs_u[idx_j, seg_j, idx_i, seg_i] = rho_si * alpha_sj
                        rho_pairs_u[idx_i, seg_i, idx_j, seg_j] = rho_si * alpha_sj
                    end
                end
            end
        end
    end
end