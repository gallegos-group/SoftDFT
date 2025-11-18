include("objective_dft.jl")


# Solution variables: rho_beads_K, rho_bonds_K, PsiC
function solver_dft(dft_system :: IsingDFT)
    
    @unpack bulk_system, geometry, fields, functionals, coulomb, numerics = dft_system

    problem_data = (bulk_system, geometry, fields, functionals, coulomb)

    objective_function = objective_dft

    solution_variables = (fields.excess.mu_ex_K, fields.excess.lng_K, fields.excess.Psi, fields.excess.PsiC)

    constraints = (identity, identity, identity, identity)

    return solver_function(
                            AndersonMixing(), 
                            problem_data,
                            objective_function,
                            solution_variables,
                            numerics,
                            constraints
                            )
end