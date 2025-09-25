# fixed_point/FixedPointSolver.jl

"""
Fixed-point solvers for nonlinear equations of the form `x = f(x)`.

This file defines the abstract interface `AbstractFixedPointSolver` and includes
specific implementations:

- `PicardMixing`: Basic adaptive fixed-point iteration
- `AndersonMixing`: Anderson acceleration using residual history and mixing

All solvers should implement a `solver_function(::SolverType, ...)` method.
"""

# Abstract interface
abstract type AbstractFixedPointSolver <: AbstractNonlinearSolver end

# Concrete solver types and their implementations
struct PicardMixing <: AbstractFixedPointSolver end
include("PicardMixing.jl")

struct AndersonMixing <: AbstractFixedPointSolver end
include("AndersonMixing.jl")
