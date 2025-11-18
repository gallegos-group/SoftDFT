using YAML
using LinearAlgebra
using LazyArrays
using IterativeSolvers

# === Solver ====
include("../NonlinearSolvers/NonlinearSolvers.jl")

# === Core Components ===
include("../molecular_system/molecular_system.jl")   # Defines `molecular_system(...)`
include("structures/IsingLST.jl")                    # Defines `iLST_Structure`
include("setup/setup_iLST.jl")                       # Builds full iLST structure
include("evaluation/eval_bulkstate.jl")               # Solves the bulk problem

# === Main Entry Point ===

"""
    iLST(input_file::String) -> IsingLST

Construct and solve an iLST system from a YAML input file.
Returns a fully initialized `IsingLST` object.
"""
function iLST(input_file::String)
    
    dataset = YAML.load_file(input_file)

    return iLST(dataset)
end

"""
    iLST(dataset::Dict) -> IsingLST

Construct and solve an iLST system from a pre-parsed YAML dataset.
"""
function iLST(dataset::Dict{Any,Any})

    bulk_system = setup_iLST(dataset)

    return iLST(bulk_system)
end

"""
    iLST(struct_iLST::IsingLST) -> IsingLST

Solve the bulk problem for a preconstructed `IsingLST`.
"""
function iLST(bulk_system::IsingLST)

    eval_bulkstate(bulk_system)

    return bulk_system
end
