
"""
    evaluate_partition_function!(
        u::Int,
        gP::Array{Float64, 3},
        gC::Array{Float64, 4},
        molsys::MolecularSystem,
        bulk::BulkState,
        species_model::freely_jointed,
    ) -> Xi

Computes the partition function `Xi` for a single chain (species `u`) in the bulk
system using the root segment (assumed to be index 1).

Arguments:
- `u`: Index of the configuration/species to evaluate.
- `gP`: Parent propagator array (size: states × states × monomers).
- `gC`: Child propagator array (size: states × states × max_children × monomers).
- `molsys`: The molecular system, including topology and properties.
- `bulk`: Thermodynamic state information (densities, chemical potentials).
- `species_model`: Chain model type (used for dispatch, currently unused).

"""
function evaluate_partition_function!(
    u::Int,
    gP::Array{Float64, 3},
    gC::Array{Float64, 4},
    molsys::MolecularSystem,
    bulk::BulkState,
    species_model::freely_jointed
)
    @unpack configurations, properties = molsys
    @unpack mu_ex = bulk
    @unpack delta_muH = properties.monomers

    config = configurations[u]
    @unpack sequence, state_family, topology = config
    @unpack children = topology

    bulk.Xi[u] = 0.0
    root = 1  # Assumes segment 1 is root

    for (idx_i, state_i) in enumerate(state_family[root])
        λ = mu_ex[state_i] + delta_muH[state_i]
        temp = exp(-λ) * sum(@views gP[:, idx_i, root])

        for (idx_c, child) in enumerate(children[root])
            temp *= sum(@views gC[:, idx_c, idx_i, root])
        end

        bulk.Xi[u] += temp
    end
end
