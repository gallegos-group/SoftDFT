"""
    solver_function(
        solver_type::AndersonMixing,
        problem_data::Tuple,
        objective_function::Function,
        solution_variables::Tuple,
        numerical_details::Dict = Dict(),
        constraints = nothing,
        error_file::String = "error_log.txt"
    ) -> Int

Solve a fixed-point problem using Anderson Mixing, optionally falling back to Picard mixing during early iterations
or when the residual is large.

### Arguments
- `solver_type::AndersonMixing`: Placeholder to indicate the solver strategy (not currently used internally).
- `problem_data::Tuple`: A tuple of inputs passed to the `objective_function`.
- `objective_function::Function`: A function that takes `problem_data` and returns a tuple of the same structure
  as `solution_variables`, representing the next guess.
- `solution_variables::Tuple`: A tuple of mutable containers (e.g., arrays) that will be updated in-place
  to hold the current solution.
- `numerical_details::Dict`: Optional dictionary specifying solver parameters:
    - `"max_iters"`: Maximum number of iterations (default = `100000`)
    - `"damping"`: Base damping factor for mixing (default = `0.5`)
    - `"mixing_max"`: Maximum allowed mixing coefficient (default = `1.0`)
    - `"switch_tole"`: Residual threshold to switch from Picard to Anderson (default = `1e-2`)
    - `"tole"`: Convergence tolerance for the residual (default = `1e-6`)
    - `"anderm"`: Depth of Anderson memory (default = `5`)
- `constraints`: A tuple of functions applied to each component of the solution vector to enforce constraints
  (e.g., non-negativity). If `nothing`, all constraints default to `identity`.
- `error_file::String`: Path to a file for logging error values at each iteration (appends to file).

### Returns
- `Int`: The number of iterations performed before convergence.

### Notes
- All entries in `solution_variables` **must be mutable** (e.g., arrays). Scalars are not allowed since they
  cannot be updated in-place inside tuples.
- This function mixes `solution_variables` toward `objective_function(problem_data)` using Picard mixing
  initially, and switches to Anderson mixing once the residual falls below `switch_tole`.
- The flattening and application of updates is handled via helper functions `flatten_solution!` and
  `apply_flat_solution!`.

### Example
```julia
solution = (zeros(10), zeros(10))
iters = solver_function(
    AndersonMixing(),
    (params,),
    my_objective,
    solution,
    Dict("tole" => 1e-8, "anderm" => 4),
)
"""

function solver_function(
    solver_type::AndersonMixing,
    problem_data::Tuple,
    objective_function::Function,
    solution_variables::Tuple,
    numerical_details::Dict = Dict(),
    constraints=nothing,
    error_file::String = "error_log.txt"
)

    ni = map(count_total_entries, solution_variables)
    n = sum(ni)

    if constraints === nothing
        constraints = fill(identity, length(solution_variables))
    end

    # --- Numerical Options ---
    max_iters   = Int(get(numerical_details, "max_iters", 100000))
    damping     = get(numerical_details, "damping", 0.5)
    mix_max     = get(numerical_details, "mixing_max", 1.0)
    switch_tole = get(numerical_details, "switch_tole", 1e-2)
    tole        = get(numerical_details, "tole", 1e-6)
    anderm      = Int(get(numerical_details, "anderm", 5))

    # --- Initialization ---
    X      = zeros(Float64, n, anderm)
    G      = zeros(Float64, n, anderm)
    X_k    = zeros(Float64, n)
    G_k    = zeros(Float64, n)
    tempX  = zeros(Float64, n, anderm - 1)
    tempG  = zeros(Float64, n, anderm - 1)
    lambdas = zeros(Float64, anderm - 1)

    err = 1.0
    iter = 0
    use_picard = true
    start_time = time()

    # open(error_file, "a") do file
    io = open(error_file, "w")
    try
        while iter < max_iters
            iter += 1
            current = objective_function(problem_data)

            col = mod(iter - 1, anderm) + 1
            flatten_solution!(X_k, solution_variables)
            flatten_solution!(G_k, current)

            @inbounds for k in eachindex(X_k)
                G_k[k] -= X_k[k]
                X[k, col] = X_k[k]
                G[k, col] = G_k[k]
            end

            err = max(maximum(abs, G_k), 1e-20)  # Prevent division by zero

            amix = min(mix_max, damping / err)

            # Log error
            msg = "Error is $err at iteration $iter.\n"
            write(io, msg)

            elapsed = time() - start_time
            print("\r$msg Elapsed: $(round(elapsed, digits=2))s"); flush(stdout)

            if iter < anderm || (err > switch_tole && use_picard)
                # Use Picard mixing
                apply_mixing!(solution_variables, current, amix)
            else
                use_picard = false

                # Build ΔX, ΔG history
                i = 0
                @inbounds for j in 1:anderm
                    if j == col
                        continue
                    end
                    i += 1
                    for k in eachindex(X_k)
                        tempX[k, i] = X[k, j] - X_k[k]
                        tempG[k, i] = G[k, j] - G_k[k]
                    end
                end

                lsqr!(lambdas, tempG, -G_k)

                result = X_k + tempX * lambdas + amix * (G_k + tempG * lambdas)

                # # Apply constraints to `result` before using it
                # idx_start = 1
                # for j in eachindex(ni)
                #     idx_end = idx_start + ni[j] - 1
                #     @. @views result[idx_start:idx_end] = constraints[j](result[idx_start:idx_end])
                #     idx_start = idx_end + 1
                # end

                # # Apply full update
                apply_mixing!(X_k, result, 1.0)
                apply_flat_solution!(solution_variables, X_k)
            end
    # println(err)
            if err < tole
                break
            end
        end
    finally
        close(io)
    end

    # println()  # Newline after progress bar
    return iter
end