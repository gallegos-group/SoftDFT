"""
    bonding_constraint(configurations::Vector, properties::PropertiesStruct)

Dispatches bonding constraints for each configuration based on its associated 
species model. Each specific species model is expected to provide its own 
method definition for:

    bonding_constraint(config::ConfigurationStruct, properties::PropertiesStruct, model::MyModel)

If no method exists for the given `model`, an error is raised.
"""
function bonding_constraint(configurations, properties)
    for (u, config) in enumerate(configurations)
        bonding_constraint(config, properties, config.species_model)
    end
end

"""
    bonding_constraint(config::ConfigurationStruct, properties::PropertiesStruct, model::AbstractSpecies)

Fallback method that raises an error when no implementation of bonding 
constraints is defined for the provided species model.
"""
function bonding_constraint(config, properties, model)
    error("No species model found for $(typeof(model)). You must define bonding_constraint(::ConfigurationStruct, ::PropertiesStruct, ::$(typeof(model)))")
end
