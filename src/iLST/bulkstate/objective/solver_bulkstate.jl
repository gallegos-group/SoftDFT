"""
    solver_bulkstate(bulk_system::IsingLST)

Solves the bulk state equations for the given `IsingLST` using Anderson mixing,
adjusting the full vector of species densities until all constraints are satisfied.

This function:
- Uses the densities in `bulk.rho.species` as the solution variables.
- Calls an objective function that enforces:
  - pSpecies chemical potential constraints
  - fixed density constraints
  - chemical equilibrium of specified reactions
  - global charge neutrality (if valences are present)

Validation:
- Calls `check_bulk_residual_consistency` to ensure the number of residuals
  matches the number of adjustable variables and that all species are properly constrained.

Applies element-wise non-negativity constraints via `x -> max(x, 0.0)`.

Modifies:
- `bulk.rho.species` in place
- Updates all dependent thermodynamic fields via the objective function

Returns:
- Solver result (if applicable), or `nothing`
"""

include("residual_bulkstate.jl")

function solver_bulkstate(bulk_system::IsingLST)

    molsys = bulk_system.molsys
    bulk   = bulk_system.bulk

    # Problem setup
    problem_data = (molsys, bulk)
    objective_function = residual_bulkstate

    # solution_variables = (bulk.rho.species, )
    solution_variables = assign_solution_variables(molsys, bulk)
    problem_data = (problem_data..., solution_variables, )

    # Merge user-specified and default solver settings
    numerics = merge(
        Dict(
            "mixing_max" => 0.1,
            "tole"       => 1e-10,
            "damping"    => 0.1
        ),
        bulk_system.numerics
    )

    # Constraints (scalar functions applied element-wise within solver_function)
    constraints = (Identity, Identity, Identity)

    return solver_function(
        AndersonResidual(),
        problem_data,
        objective_function,
        solution_variables,
        numerics,
        constraints
    )
end


function assign_solution_variables(molsys, bulk)
    @unpack species, input_pSpecies, adjustable = molsys.properties.species
    @unpack rho = bulk
    
    reactions = molsys.properties.reactions

    Q_solver = [0.0]  # scalar as a 1-element array for mutability

    log_rho_rxn = Float64[]
    log_rho_pH  = Float64[]

    # Identify relevant species from real reactions (exclude bookkeeping)
    species_in_reactions = Set{String}()
    for rxn in reactions
        if all(c == 0 for c in rxn.coeffs)
            continue
        end
        for s in rxn.species
            push!(species_in_reactions, s)
        end
    end

    for (i, s) in enumerate(species)
        ρi = rho.species[i]

        if !isnan(input_pSpecies[i])
            # Species constrained via pH
            push!(log_rho_pH, log(ρi))

        elseif s in species_in_reactions
            # Species appears in a reaction
            push!(log_rho_rxn, log(ρi))
        else
            # Fixed and not used in anything relevant → do not include
            continue
        end
    end

    return (log_rho_rxn, log_rho_pH, Q_solver, )
end