"""
    solver_function(
        solver_type::AndersonResidual,
        problem_data::Tuple,
        residual_function::Function,
        solution_variables::Tuple,
        numerical_details::Dict = Dict(),
        constraints = nothing,
        error_file::String = "error_log.txt"
    ) -> Int

Solve a nonlinear system `f(x) = 0` using Anderson mixing on the residual.

### Arguments
- `solver_type::AndersonResidual`: Dispatch tag for Anderson mixing applied to residuals.
- `problem_data::Tuple`: Static arguments passed to the residual function.
- `residual_function::Function`: Function of form `residual = f(problem_data...)` that returns a residual with same structure as `solution_variables`.
- `solution_variables::Tuple`: Mutable containers that will be updated in-place.
- `numerical_details::Dict`: Solver parameters:
    - `"max_iters"`, `"damping"`, `"mixing_max"`, `"tole"`, `"anderm"` (same as fixed-point version)
- `constraints`: Functions applied to each element of the solution vector (default = identity).
- `error_file`: Log file for error per iteration.

### Returns
- `Int`: Number of iterations used
"""
function solver_function(
    solver_type::AndersonResidual,
    problem_data::Tuple,
    residual_function::Function,
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

    max_iters   = get(numerical_details, "max_iters", 100000)
    damping     = get(numerical_details, "damping", 0.1)
    mix_max     = get(numerical_details, "mixing_max", 0.1)
    switch_mix  = get(numerical_details, "switch_mix", 0.01)
    tole        = get(numerical_details, "tole", 1e-6)
    anderm      = get(numerical_details, "anderm", 5)

    # Storage for Anderson mixing
    X      = zeros(Float64, n, anderm)
    F      = zeros(Float64, n, anderm)
    X_k    = zeros(Float64, n)
    F_k    = zeros(Float64, n)
    tempX  = zeros(Float64, n, anderm - 1)
    tempF  = zeros(Float64, n, anderm - 1)
    lambdas = zeros(Float64, anderm - 1)

    # Initialize Anderson memory with initial guess
    flatten_solution!(X_k, solution_variables)
    flatten_solution!(F_k, residual_function(problem_data))

    @. X[:, 1] = X_k
    @. F[:, 1] = F_k

    iter = 0
    err = 1.0
    start_time = time()
        open(error_file, "w") do file
    while iter < max_iters
        iter += 1

        # Apply current guess to solution_variables
        apply_flat_solution!(solution_variables, X_k)

        # Residual function
        residual = residual_function(problem_data)

        # Flatten residual
        flatten_solution!(F_k, residual)

        # Compute error and mixing
        err = max(maximum(abs, F_k), 1e-20)
        amix = min(mix_max, damping / err)

        # Store history: guess → residual
        col = mod(iter, anderm) + 1
        @. @views X[:, col] = X_k
        @. @views F[:, col] = F_k

        # Logging
        msg = "Residual error is $err at iteration $iter.\n"

            write(file, msg)

        elapsed = time() - start_time
        print("\r$msg Elapsed: $(round(elapsed, digits=2))s"); flush(stdout)

        if iter < anderm || err > switch_mix
            prev_col = mod(iter-1, anderm) + 1

            ΔX = X_k .- X[:, prev_col]
            ΔF = F_k .- F[:, prev_col]

            # Project F_k onto ΔF
            numer = dot(F_k, ΔF)
            denom = dot(ΔF, ΔF)
            λ = denom == 0.0 ? 0.0 : numer / denom

            # Picard update with projection and damping
            @. X_k = X_k - λ * ΔX - amix * (F_k - λ * ΔF)
        else
            i = 0
            for j in 1:anderm
                if j == col
                    continue
                end
                i += 1
                @. @views tempX[:, i] = X[:, j] - X_k
                @. @views tempF[:, i] = F[:, j] - F_k
            end

            lsqr!(lambdas, tempF, -F_k)

            result = X_k + tempX * lambdas + amix*(F_k + tempF*lambdas)

            @. X_k = result
        end

        if err < tole
            break
        end
    end

    end
    # println()
    return iter
end

