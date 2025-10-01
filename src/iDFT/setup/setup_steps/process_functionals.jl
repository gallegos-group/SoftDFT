include("../../functionals/Abstract_Types/functional_registry.jl")

"""
    process_functionals(models_list::Vector{Symbol}, bulk_system::IsingLST, geometry::CoordSystem)

Given a list of functional names, this function:
- Validates the names using the model registry,
- Retrieves and constructs each functional,
- Returns a list of initialized `AbstractFunctional` objects.

Arguments:
- `models_list` — list of functional names as Symbols (e.g., `[:ideal, :mf]`)
- `bulk_system` — contains `molsys` and `bulk`, used for initialization
- `geometry` — the grid/coordinate system

Returns:
- `Vector{AbstractFunctional}`
"""
function process_functionals(bulk_system::IsingLST, geometry::CoordSystem)
    models_list = bulk_system.molsys.properties.fe_model

    validate_functional_keys(models_list)

    functionals = Vector{AbstractFunctional}()
    for name in models_list
        functional_type = get_functional(name)
    
        if !is_ideal(functional_type)
            model_instance = construct_functional(functional_type, bulk_system.molsys, bulk_system.bulk, geometry)
            push!(functionals, model_instance)
        end
    end

    return functionals
end

