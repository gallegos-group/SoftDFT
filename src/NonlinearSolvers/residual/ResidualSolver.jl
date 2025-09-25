# residual/ResidualSolver.jl

"""
Residual-based solvers for systems of the form `f(x) = 0`.

Defines:
- `AbstractResidualSolver`: base type for residual solvers
- `PicardResidual`: simple residual-based descent update
- `AndersonResidual`: accelerated residual-based Anderson mixing

Each solver dispatches to a method `solver_function(::SolverType, ...)`.
"""

# Abstract interface
abstract type AbstractResidualSolver <: AbstractNonlinearSolver end

# Concrete solver types
struct PicardResidual <: AbstractResidualSolver end
include("PicardResidual.jl")

struct AndersonResidual <: AbstractResidualSolver end
include("AndersonResidual.jl")
