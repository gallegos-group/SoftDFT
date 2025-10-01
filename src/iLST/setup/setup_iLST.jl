# =========================================================================
# setup_iLST.jl — Constructs IsingLST for liquid-state theory models
# =========================================================================


include("setup_steps/parse_bulkstate.jl")
include("setup_steps/parse_numerical.jl")
include("setup_steps/process_evaluation.jl")

# -------------------------------------------------------------------------
# Top-level iLST setup function
# -------------------------------------------------------------------------

"""
    setup_iLST(data::Dict, molsys::MolecularSystem) -> IsingLST

Initializes an `IsingLST` used in self-consistent field or liquid-state
theory calculations.

Steps:
1. Processes densities and bond topology.
2. Constructs model components from the model list in `data["model"]`.
3. Builds the `BulkState` (`bulk`) with initialized thermodynamic fields.
4. Extracts or defaults numerical solver parameters.

Returns:
- A fully constructed `IsingLST` ready for in-place iterative solution.
"""

function setup_iLST(dataset::Dict{Any, Any}, molsys::MolecularSystem)

    evaluation = process_evaluation(dataset)

    bulk = parse_bulkstate(molsys, evaluation)

    numerics = parse_numerical(dataset)

    return IsingLST(molsys, bulk, numerics)
end


function setup_iLST(dataset::Dict{Any, Any})

    molsys = molecular_system(dataset)

    evaluation = process_evaluation(dataset)

    bulk = parse_bulkstate(molsys, evaluation)
      
    numerics = parse_numerical(dataset)

    process_evaluation(dataset)

    return IsingLST(molsys, bulk, numerics)
end