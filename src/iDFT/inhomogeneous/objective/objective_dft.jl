include("functions/update_spatial_structures.jl")
include("functions/eval_spatial_density.jl")
include("functions/eval_electroneutrality.jl")
include("functions/eval_rho_hat.jl")
include("functions/eval_mu_functional.jl")

function objective_dft(dft_system :: IsingDFT)
    @unpack bulk_system, geometry, fields, functionals = dft_system

    return objective_dft(bulk_system, geometry, fields, functionals)
end

function objective_dft(problem_data::Tuple)
    bulk_system = get_first_of_type(IsingLST, problem_data)
    geometry    = get_first_of_type(CoordSystem, problem_data)
    fields      = get_first_of_type(SpatialFields, problem_data)
    functionals = get_first_of_type(Vector{<:AbstractFunctional}, problem_data)

    return objective_dft(bulk_system, geometry, fields, functionals)
end


function objective_dft(bulk_system, geometry, fields, functionals)

    update_spatial_structures!(bulk_system, geometry, fields)

    eval_rho_hat!(bulk_system, geometry, fields)

    eval_mu_functional!(bulk_system, geometry, fields, functionals)

    rho1_segments_K, rho1_bonds_K = 
        eval_spatial_density(bulk_system, geometry, fields)

    Psi1C = eval_electroneutrality(bulk_system, geometry, fields)

    return rho1_segments_K, rho1_bonds_K, Psi1C
end
