

include("initialize_densities.jl")
include("initialize_fe_model.jl")



"""
    parse_bulkstate(dataset, molsys) -> BulkState

Builds a `BulkState` object from the molecular system and model list
specified in `dataset`.

Steps:
1. Computes initial densities using `process_densities`.
2. Constructs model components from 'process_model'.
3. Returns a fully initialized `BulkState` ready for iterative solution.
"""
function parse_bulkstate(molsys :: MolecularSystem, evaluation :: F) where {F <: Tuple{Vararg{AbstractEvaluation}}}
    rho = initalize_densities(molsys)
    fe_model = initialize_fe_model(molsys.properties.fe_model, rho)

    return BulkState(rho, fe_model, evaluation)
end