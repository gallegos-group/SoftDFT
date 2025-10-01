

"""
    compute_parent_propagators(
        u::Int,
        gC::Array{Float64, 4},
        molsys::MolecularSystem,
        bulk::BulkState,
        species_model::freely_jointed
    ) -> gP

Computes the parent (backward) propagators `gP` for a structured chain (species `u`)
using precomputed child propagators `gC` and the current thermodynamic state.

Arguments:
- `u`: Index of the species/configuration.
- `gC`: Child propagator array (size: states × states × max_children × monomers).
- `molsys`: The molecular system definition, including configurations and monomer properties.
- `bulk`: Thermodynamic quantities including chemical potentials and bond topology.
- `species_model`: Chain model type (currently used for dispatch only).

Returns:
- `gP`: The computed parent propagator array (size: states × states × monomers).
"""
function compute_parent_propagators(
    u::Int,
    gC::Array{Float64, 4},
    molsys::MolecularSystem,
    bulk::BulkState,
    species_model::freely_jointed
)

    @unpack configurations, properties = molsys
    @unpack sequence, state_family, topology = configurations[u]
    @unpack parents, children, levels = topology
    @unpack mu_ex, lng = bulk
    @unpack delta_muH = properties.monomers
    @unpack bond_types = properties

    MTS = length(sequence)
    max_states = max(maximum(x -> length(x), state_family), 1)
    gP = zeros(Float64, max_states, max_states, MTS)

    for level in eachindex(levels)
        for child in levels[level]
            v = sequence[child]

            for (idx_j, state_j) in enumerate(state_family[child])
                gP[1, idx_j, child] = 1.0  # base case

                for parent in parents[child]
                    for (idx_i, state_i) in enumerate(state_family[parent])
                        pair = state_i <= state_j ? (state_i, state_j) : (state_j, state_i)
                        bond = findfirst(x -> x == pair, bond_types)
                        
                        val = exp(-delta_muH[state_i] - mu_ex[state_i] + lng[bond])
                        val *= sum(@views gP[:, idx_i, parent])

                        for (idx_oc, other_child) in enumerate(children[parent])
                            if other_child != child
                                val *= sum(@views gC[:, idx_oc, idx_i, parent])
                            end
                        end

                        gP[idx_i, idx_j, child] = val
                    end
                end
            end
        end
    end

    return gP
end
