using YAML
using FFTW
using LinearAlgebra
using LazyArrays
using LoopVectorization
using IterativeSolvers
using IterTools

# === Upstream Dependencies ===
include("../iLST/iLST.jl")  # Provides IsingLST: `molsys`, `bulk_state`
include("../geometry/geometry.jl") # Provides 'geometry'

# === Core Components ===
include("structures/IsingDFT.jl")
include("setup/setup_iDFT.jl")
include("evaluation/eval_inhomogeneous.jl")

"""
    iDFT(input_file::String) -> IsingDFT

Construct and solve an inhomogeneous DFT system from a YAML input file.
"""
function iDFT(input_file::String)
    dataset = YAML.load_file(input_file)

    # Solve bulk state (includes mol_sys and bulk_state)
    bulk_system = iLST(dataset)

    # Build geometry
    geom = geometry_system(dataset)

    # Build geometry, fields, functionals
    dft_system = setup_iDFT(dataset, bulk_system, geom)

    # Solve inhomogeneous problem
    eval_inhomogeneous(dft_system)

    return dft_system
end

function iDFT(input_file::String, a)
    dataset = YAML.load_file(input_file)

    # Solve bulk state (includes mol_sys and bulk_state)
    bulk_system = iLST(dataset)

    # Build geometry
    geom = geometry_system(dataset)

    # Build geometry, fields, functionals
    dft_system = setup_iDFT(dataset, bulk_system, geom)

    # Solve inhomogeneous problem
    eval_inhomogeneous(dft_system, a)

    return dft_system
end