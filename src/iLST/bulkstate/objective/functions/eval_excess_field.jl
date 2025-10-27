"""
    excess_field!(molsys, bulk)

Computes the total chemical potential, contact values (lng), and
species-level chemical potentials by looping over all free energy models.

Clears all fields first using `@.` assignment, then accumulates contributions.
"""
function eval_excess_field!(molsys :: MolecularSystem, bulk :: BulkState)

    @unpack mu_ex, lng = bulk

    # Zero out excess fields
    @. mu_ex       = 0.0
    @. lng         = 0.0

    # Accumulate contributions from each model term
    # Explicitly unroll for known model count
    ntuple(length(bulk.fe_model)) do i
        model_term = bulk.fe_model[i]
        if !is_ideal(model_term)
            chemical_potential!(molsys, bulk, model_term)
        end
    end

    return nothing
end