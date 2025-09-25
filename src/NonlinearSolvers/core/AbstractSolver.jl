# core/AbstractSolver.jl

abstract type AbstractNonlinearSolver end

include("SolverUtils.jl")

"""
    solve(solver::AbstractNonlinearSolver, f, x0; kwargs...) → x

Main entry point for solving nonlinear problems.
- If `solver` is a fixed-point solver: solves `x = f(x)`
- If `solver` is a residual solver: solves `f(x) = 0`
"""
function solve(solver::AbstractNonlinearSolver, f, x0; kwargs...)
    error("solve not implemented for solver of type $(typeof(solver))")
end