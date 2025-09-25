include("objective/solver_dft.jl")
include("analysis/analysis.jl")

function eval_inhomogeneous(dft_system :: IsingDFT)

    eval_external_field(dft_system)

    process_guess(dft_system)

    solver_dft(dft_system)

    call_analysis(dft_system)

    println("inhomogeneous system solved.\n")
end

function eval_inhomogeneous(dft_system :: IsingDFT, a :: Int)

    eval_external_field(dft_system)

    # initial_guess(dft_system.bulk_system, dft_system.geometry, dft_system.fields)

    solver_dft(dft_system)

    call_analysis(dft_system)

    println("inhomogeneous system solved.\n")
end

