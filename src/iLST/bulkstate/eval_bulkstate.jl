include("objective/solver_bulkstate.jl")
include("analysis/analysis.jl")

"""
    eval_bulkstate(struct_iLST::IsingLST)

Solves the bulk self-consistency equations and updates all internal fields in `bulk`.

Steps:
1. Solves model-specific bulk equations.
2. Computes and updates the bulk pressure.

Returns:
- Nothing (all updates are in-place in `struct_iLST.bulk`).
"""
function eval_bulkstate(bulk_system :: IsingLST)
    solver_bulkstate(bulk_system)
    
    eval_bulk_pressure(bulk_system.molsys, bulk_system.bulk)

    # println("Bulk system solved.\n")
end