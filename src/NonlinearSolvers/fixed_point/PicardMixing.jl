"""
    solver_function(
        solver_type::AdaptivePicardMixing,
        problem_data::Tuple,
        objective_function::Function,
        solution_variables::Tuple,
        numerical_details::Dict = Dict(),
        constraints = nothing,
        error_file::String = "error_log.txt"
    ) -> Tuple

Solve a fixed-point problem using adaptive Picard mixing.

At each iteration, the solution is updated toward the output of `objective_function` using
an adaptive damping factor based on the current residual. Convergence is declared when
the residual norm falls below the specified tolerance.

### Arguments
- `solver_type::AdaptivePicardMixing`: A type tag used to dispatch the appropriate solver logic.
- `problem_data::Tuple`: Static input data passed to `objective_function` each iteration.
- `objective_function::Function`: A function that accepts `problem_data` and returns a tuple of values
  matching the structure of `solution_variables`.
- `solution_variables::Tuple`: A tuple of mutable objects (e.g. arrays) that are updated in-place.
- `numerical_details::Dict`: Optional dictionary of solver parameters:
    - `"max_iters"`: Maximum number of iterations (default = `1_000_000`)
    - `"damping"`: Base damping factor (default = `0.5`)
    - `"mixing_max"`: Maximum allowed damping (default = `1.0`)
    - `"tole"`: Residual tolerance for convergence (default = `1e-6`)
- `constraints`: Tuple of constraint functions applied to each element of `solution_variables`.
  If not provided, defaults to identity functions.
- `error_file::String`: Path to file for logging the error at each iteration (appends if file exists).

### Returns
- `Tuple`: The final `solution_variables`, updated in-place and also returned.

### Notes
- All entries in `solution_variables` must be mutable (e.g., arrays). Scalars are **not** supported.
- Mixing uses `amix = min(mixing_max, damping / err)` to adaptively scale updates.
- If convergence is not reached within `max_iters`, an error is thrown.
- Residual is computed as the maximum absolute difference between current and previous solution values.

### Example
```julia
solution = (zeros(100), zeros(100))
iters = solver_function(
    AdaptivePicardMixing(),
    (params,),
    my_objective,
    solution,
    Dict("tole" => 1e-8, "damping" => 0.3)
)
"""

function solver_function(
    solver_type::PicardMixing,
    problem_data::Tuple,
    objective_function::Function,
    solution_variables::Tuple,
    numerical_details::Dict = Dict(),
    constraints=nothing,
    error_file::String = "error_log.txt"
)

    ni = length.(solution_variables)
    n = sum(ni)

    # Initialize constraints to identity if not given
    if constraints === nothing
        constraints = fill(identity, length(solution_variables))
    end

    # Extract numerical options
    max_iters = get(numerical_details, "max_iters", 1000000)
    damping   = get(numerical_details, "damping", 0.5)
    mix_max   = get(numerical_details, "mixing_max", 1.0)
    tole      = get(numerical_details, "tole", 1e-6)

    iter = 0
    err = 1.0
    current = nothing

    ni = map(count_total_entries, solution_variables)
    n = sum(ni)

    X_k    = zeros(Float64, n)
    G_k    = zeros(Float64, n)


    while iter < max_iters && err > tole
        iter += 1
        current = objective_function(problem_data)

        # Compute residual error
        # err = 0.0
        # for i in eachindex(solution_variables)
        #     x = solution_variables[i]
        #     y = current[i]
        #     @inbounds for idx in eachindex(x)
        #         err = max(err, abs(y[idx] - x[idx]))
        #     end
        # end

        flatten_solution!(X_k, solution_variables)
        flatten_solution!(G_k, current)
        err = max(maximum(abs, G_k), 1e-20)

        # Adaptive mixing coefficient
        amix = min(mix_max, damping / err)
        apply_mixing!(solution_variables, current, amix)
        
        # Apply mixing with constraints
        # for i in eachindex(solution_variables)
        #     if solution_variables[i] isa AbstractArray
        #         @. solution_variables[i] = constraints[i]((1.0 - amix) * solution_variables[i] + amix * current[i])
        #     else
        #         error("Not an abstract array")
        #     end
        # end

        # Logging
        msg = "Error is $err at iteration $iter.\n"
        println(msg)
        open(error_file, "a") do file
            write(file, msg)
        end
    end

    if err > tole
        error("Solver failed to converge within $max_iters iterations. Final error = $err")
    end

    return solution_variables
end