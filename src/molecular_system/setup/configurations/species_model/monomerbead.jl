"""
    monomerbead <: AbstractSpecies

Represents a standalone, non-polymeric species with no chain connectivity.
Used for modeling individual molecules without internal bonding.
"""
struct monomerbead <: AbstractSpecies end


"""
    bonding_constraint(config::ConfigurationStruct, properties::PropertiesStruct, chain_model::monomerbead)

No bonding constraints are imposed for the `monomerbead` model.
This function is a no-op, included for compatibility with the constraint dispatch system.
"""
function bonding_constraint(config, properties, chain_model :: monomerbead)
    # intentionally does nothing
end
