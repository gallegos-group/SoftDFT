include("../../free_energy/Abstract_Types/model_registry.jl")

"""
    initalize_fe_model(models_list::Vector{Symbol}, rho)

Given a list of model names as symbols and a density structure `rho`, this function:
1. Validates the provided model names against the `ModelRegistry`,
2. Retrieves the constructor for each model,
3. Instantiates each model using `rho`,
4. Returns a vector of model structures representing the full thermodynamic model.

Arguments:
- `models_list::Vector{Symbol}` — names of models to include (e.g., `[:ideal, :sw]`)
- `rho` — density structure used to initialize the models

Returns:
- `Vector{AbstractFreeEnergy}` — list of model components

Throws:
- An error if any model name is not found in the registry.
"""
function initialize_fe_model(models_list::Vector{Symbol}, rho)
    validate_model_keys(models_list)

    model = Vector{AbstractFreeEnergy}()
    for name in models_list
        model_type = get_model(name)
        model_instance = construct_free_energy(model_type, rho)
        push!(model, model_instance)
    end

    return model
end

