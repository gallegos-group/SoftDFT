include("objective/solver_dft.jl")
include("analysis/analysis.jl")

function eval_inhomogeneous(dft_system :: IsingDFT)

    solver_dft(dft_system)

    call_analysis(dft_system)

    println("\ninhomogeneous system solved.\n")
end
