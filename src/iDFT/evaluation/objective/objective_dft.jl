include("functions/update_spatial_structures.jl")
include("functions/eval_spatial_density.jl")
include("functions/eval_rho_hat.jl")
include("functions/eval_mu_functional.jl")

function objective_dft(problem_data::Tuple)
    bulk_system = get_first_of_type(IsingLST, problem_data)
    geometry    = get_first_of_type(CoordSystem, problem_data)
    fields      = get_first_of_type(SpatialFields, problem_data)
    functionals = get_first_of_type(Tuple, problem_data)
    coulomb     = get_first_of_type(AbstractCoulomb, problem_data)

    return objective_dft(bulk_system, geometry, fields, functionals, coulomb)
end


function objective_dft(bulk_system, geometry, fields, functionals, coulomb)

    eval_spatial_density!(bulk_system, geometry, fields)
    
    update_spatial_structures!(bulk_system, geometry, fields)

    eval_rho_hat!(bulk_system, geometry, fields)

    mu_ex_K, lng_K = eval_mu_functionals(bulk_system, geometry, fields, functionals)

    Psi, PsiC = solve_coulomb(bulk_system, geometry, fields, coulomb)

    return mu_ex_K, lng_K, Psi, PsiC
end
