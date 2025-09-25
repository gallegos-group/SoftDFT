"""
    bulk_chemical_potential(molsys, bulk)

Computes the total chemical potential, contact values (lng), and
species-level chemical potentials by looping over all free energy models.

Clears all fields first using `@.` assignment, then accumulates contributions.
"""
function bulk_chemical_potential!(molsys :: MolecularSystem, bulk :: BulkState)

    @unpack mu_species, mu_ex, lng = bulk

    # Zero out chemical potential fields
    @. mu_species  = 0.0
    @. mu_ex       = 0.0
    @. lng         = 0.0

    # Update ideal and excess fields
    for model in bulk.fe_model
        chemical_potential!(molsys, bulk, model)
    end

      # Evaluate excess contribution for each species
      for (u, config) in enumerate(molsys.configurations)
        species_model = config.species_model
        bulk_chemical_potential!(u, molsys, bulk, species_model)
    end  
end

function bulk_chemical_potential!(u, molsys :: MolecularSystem, bulk :: BulkState, species_model)
    error("Species $u has no defined species_model: $species_model.")
end
