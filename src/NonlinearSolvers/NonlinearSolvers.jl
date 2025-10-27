"""
    NonlinearSolvers.jl

Modular framework for solving nonlinear equations of the form:

- Fixed-point problems: `x = f(x)`
- Residual-based problems: `f(x) = 0`
- (Optionally) continuation problems: parameterized solutions along a curve

### Structure

This file organizes the solver types and interfaces into three main categories:

- `core/AbstractSolver.jl`: Defines the abstract base types and fallback `solve` interface.
- `fixed_point/FixedPointSolver.jl`: Registers `AbstractFixedPointSolver` and optional `solve` overload.
- `residual/ResidualSolver.jl`: Registers `AbstractResidualSolver` and optional `solve` overload.
- `continuation/Continuation.jl`: (Optional) Defines continuation algorithms for tracing solution branches.

Each solver implementation (e.g., `AndersonMixing`, `PicardResidual`) dispatches on a `solver_function(::SolverType, ...)` method that performs the actual iteration logic.

### Usage

To use a solver, define a residual or fixed-point function and pass it along with an initial guess:

```julia
using NonlinearSolvers

solution = (zeros(100),)
iters = solver_function(
    AndersonResidual(),
    (params,),
    my_residual_function,
    solution,
    Dict("tole" => 1e-8, "anderm" => 5),
)
Alternatively, you may use the generic interface:

solve(AndersonResidual(), my_residual_function, solution; numerical_details = Dict(...))

"""

using UnPack

# Abstract types and utilities
include("core/AbstractSolver.jl")

# Fixed-point solvers
include("fixed_point/FixedPointSolver.jl")

# Residual solvers
include("residual/ResidualSolver.jl")

# Continuation methods
include("continuation/Continuation.jl")

