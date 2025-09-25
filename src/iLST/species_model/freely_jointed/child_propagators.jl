"""
    compute_child_propagators(u, molsys, bulk, freely_jointed)

Compute the child propagator array `gC` for species `u` based on a freely-jointed chain model.

This function computes the forward (child) propagators recursively from the leaves up
using the segment topology and bond connectivity.

Arguments:
- `u`: Index of the species being evaluated.
- `molsys`: MolecularSystem containing sequence and topology information.
- `bulk`: BulkState containing chemical potentials and bond contributions.
- `freely_jointed`: Dispatch tag indicating chain model used.

Returns:
- `gC::Array{Float64,4}` with dimensions `(child_state, parent_state, child_index, parent_segment)`
"""
function compute_child_propagators(
    u::Int,
    molsys::MolecularSystem,
    bulk::BulkState,
    species_model::freely_jointed
)

    @unpack configurations, properties = molsys
    @unpack sequence, state_family, topology = configurations[u]
    @unpack children, levels = topology
    @unpack mu_ex, lng = bulk
    @unpack delta_muH = properties.monomers
    @unpack bond_types = properties

    MTS = length(sequence)
    max_states = max(maximum(x -> length(x), state_family), 1)
    max_childs = max(maximum(x -> length(x), children), 1)

    gC = zeros(Float64, max_states, max_childs, max_states, MTS)

    for level in reverse(eachindex(levels))  # # bottom-up: start from leaves
        for parent in levels[level]
            for (idx_i, state_i) in enumerate(state_family[parent])

                gC[1, 1, idx_i, parent] = 1.0  # base case

                for (idx_c, child) in enumerate(children[parent])
                    for (idx_j, state_j) in enumerate(state_family[child])
                        pair = state_i <= state_j ? (state_i, state_j) : (state_j, state_i)
                        bond = findfirst(x -> x == pair, bond_types)

                        gC[idx_j, idx_c, idx_i, parent] =
                            exp(-delta_muH[state_j] - mu_ex[state_j] + lng[bond])

                        for (idx_gc, grandchild) in enumerate(children[child])
                            gC[idx_j, idx_c, idx_i, parent] *=
                                sum(@views gC[:, idx_gc, idx_j, child])
                        end
                    end
                end
            end
        end
    end

    return gC
end