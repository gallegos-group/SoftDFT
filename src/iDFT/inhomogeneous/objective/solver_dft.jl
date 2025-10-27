include("objective_dft.jl")


# Solution variables: rho_beads_K, rho_bonds_K, PsiC
function solver_dft(dft_system :: IsingDFT)
    
    @unpack bulk_system, geometry, fields, functionals, numerics = dft_system

    problem_data = (bulk_system, geometry, fields, functionals, )

    objective_function = objective_dft

    solution_variables = (fields.rho_K.segments, fields.rho_K.bonds, fields.excess.PsiC)

    constraints = (identity, identity, identity)

    return solver_function(
                            AndersonMixing(), 
                            problem_data,
                            objective_function,
                            solution_variables,
                            numerics,
                            constraints
                            )
end