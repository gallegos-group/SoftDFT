
"""
    eval_residuals!(molsys::MolecularSystem, bulk::BulkState) -> (residuals,)

Evaluates the residual vector for the bulk system, enforcing the following constraints:

1. **pSpecies constraints**: For species with a defined protonation potential (`pSpecies`),
   residuals are defined as `μ_i - pK_i * log(10)`.

2. **Fixed density constraints**: For species that are not adjustable and do not have a `pSpecies` constraint,
   residuals enforce the density to remain equal to the specified `input_densities`.

3. **Reaction equilibrium constraints**: For each chemical reaction (excluding bookkeeping),
   a residual is added to enforce `∑ ν_i μ_i - pK = 0`.

4. **Charge neutrality constraint**: The total charge `∑ z_i ρ_i` is added as a final residual term.
   If no species are adjustable and the system is not neutral, an error is raised.

Arguments:
- `molsys`: Molecular system definition including species, pSpecies, and reactions
- `bulk`: Current bulk state including densities and chemical potentials

Returns:
- A 1-tuple containing the residual vector (as required by solver interface)
"""
function eval_residuals!(molsys::MolecularSystem, bulk::BulkState)
    @unpack properties = molsys
    @unpack mu_species, rho = bulk
    @unpack species, input_pSpecies, species_charge = properties.species
    reactions = properties.reactions

    residuals = Float64[]
    unit_conv = concentration_to_density(1.0, molsys.constants)

    # === 1. pSpecies-constrained residuals ===
    for (i, _) in enumerate(species)
        if !isnan(input_pSpecies[i])
            residual = mu_species[i] + input_pSpecies[i] * log(10) - log(unit_conv)
            push!(residuals, residual)
        end
    end

    # === 2. Reaction equilibrium residuals ===
    for rxn in reactions
        if all(c == 0 for c in rxn.coeffs)
            continue  # bookkeeping-only reaction
        end
        
        μsum = sum(ν * mu_species[findfirst(isequal(s), species)]
                   for (s, ν) in zip(rxn.species, rxn.coeffs))
        residual = μsum - rxn.pK[1]
        push!(residuals, residual)
    end

    # === 3. Charge neutrality residual ===
    total_charge = -sum(@. @~ species_charge * rho.species)
    push!(residuals, total_charge)
    
    return (residuals,)
end