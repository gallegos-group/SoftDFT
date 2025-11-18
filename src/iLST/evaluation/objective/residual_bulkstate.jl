include("functions/balance_charge.jl")
include("functions/update_bulk_structures.jl")
include("functions/eval_bulk_density.jl")
include("functions/eval_residuals.jl")

"""
    residual_bulkstate(problem_data::Tuple)

Extracts `BulkState` and `MolecularSystem` from the input tuple by type,
and computes the objective residuals for the bulk system.

This version supports generic solver interfaces that pass problem data as a tuple.

Returns:
- Residual vector tuple from `eval_residuals!`
"""
function residual_bulkstate(problem_data::Tuple)

    molsys = get_first_of_type(MolecularSystem, problem_data)
    bulk   = get_first_of_type(BulkState, problem_data)
    soln   = get_first_of_type(Tuple, problem_data)

    return residual_bulkstate(molsys, bulk, soln)
end

function residual_bulkstate(molsys :: MolecularSystem, bulk :: BulkState, solution)

    update_bulk_structures!(molsys, bulk, solution)

    eval_bulk_density!(molsys, bulk)

    bulk_chemical_potential!(molsys, bulk)

    return eval_residuals!(molsys, bulk)
end


function get_first_of_type(::Type{T}, tup::Tuple) where T
    idx = findfirst(x -> x isa T, tup)
    isnothing(idx) && error("Missing item of type $(T)")
    return tup[idx]
end