"""
    solver_function(
        solver_type::PicardResidual,
        problem_data::Tuple,
        residual_function::Function,
        solution_variables::Tuple,
        numerical_details::Dict = Dict(),
        constraints = nothing,
        error_file::String = "error_log.txt"
    ) -> Tuple

Solve a residual-based system `f(x) = 0` using a structure-agnostic Picard-style update.

This method is a single-step residual solver that does not assume a one-to-one correspondence
between residual and solution components. Instead of performing per-variable updates like standard
Picard iteration, it applies a projected update based on the change in residual and solution between
iterations. This makes it suitable for globally coupled systems, such as those involving electroneutrality,
mass balance, or other nonlocal constraints.

### Arguments
- `solver_type::PicardResidual`: Type tag for dispatching this solver.
- `problem_data::Tuple`: Static input arguments passed to `residual_function`.
- `residual_function::Function`: A function that returns a structure matching `solution_variables`
  representing the residual `f(x)` at each iteration.
- `solution_variables::Tuple`: A tuple of mutable containers (e.g., arrays) representing the current guess.
  These are updated in-place each iteration.
- `numerical_details::Dict`: Optional dictionary with the following keys:
    - `"max_iters"`: Maximum number of iterations (default = `100_000`)
    - `"damping"`: Base damping factor (default = `0.5`)
    - `"mixing_max"`: Maximum damping applied (default = `1.0`)
    - `"tole"`: Convergence tolerance on residual norm (default = `1e-6`)
- `constraints`: Tuple of constraint functions (same length as `solution_variables`) applied to
  slices of the flattened solution. Defaults to identity functions.
- `error_file::String`: Path to log residual error at each iteration.

### Returns
- `Tuple`: The final in-place updated `solution_variables`.

### Notes
- This method performs a projected update using the previous step’s residual and solution change.
- It is a generalization of Picard iteration that does **not** require the residual to align
  structurally with the solution vector.
- The first iteration is used to evaluate direction and does not perform an update.
- Suitable for residuals involving coupled constraints or aggregate conditions.

### Example
```julia
solution = (zeros(100),)
iters = solver_function(
    PicardResidual(),
    (params,),
    compute_residual,
    solution,
    Dict("tole" => 1e-8, "damping" => 0.2)
)
"""

function solver_function(
    solver_type::PicardResidual,
    problem_data::Tuple,
    residual_function::Function,
    solution_variables::Tuple,
    numerical_details::Dict = Dict(),
    constraints = nothing,
    error_file::String = "error_log.txt"
)

    ni = map(count_total_entries, solution_variables)
    n = sum(ni)

    if constraints === nothing
        constraints = fill(identity, length(solution_variables))
    end

    max_iters = get(numerical_details, "max_iters", 100_000)
    damping   = get(numerical_details, "damping", 0.1)
    mix_max   = get(numerical_details, "mixing_max", 0.1)
    tole      = get(numerical_details, "tole", 1e-6)

    X_prev = zeros(Float64, n)
    F_prev = zeros(Float64, n)
    X_k    = zeros(Float64, n)
    F_k    = zeros(Float64, n)

    # Initial residual
    flatten_solution!(X_k, solution_variables)
    flatten_solution!(X_prev, solution_variables)
    flatten_solution!(F_prev, residual_function(problem_data))

    err = 1.0
    iter = 1
    start_time = time()

    while iter <= max_iters && err > tole

        # Save current guess before applying (used for ΔX)
        X_k_original = copy(X_k)

        # Apply guess to internal structures
        apply_flat_solution!(solution_variables, X_k)

        # Evaluate residual (may update solution_variables)
        residual = residual_function(problem_data)

        # Flatten residual
        flatten_solution!(F_k, residual)

        # Compute error and mixing
        err = max(maximum(abs, F_k), 1e-20)
        amix = min(mix_max, damping / err)

        # Compute ΔX, ΔF
        ΔX = X_k_original .- X_prev
        ΔF = F_k .- F_prev

        # Project F_k onto ΔF
        numer = dot(F_k, ΔF)
        denom = dot(ΔF, ΔF)
        λ = denom == 0.0 ? 0.0 : numer / denom

        # Picard update with projection and damping
        @. X_k = X_k_original - λ * ΔX - amix * (F_k - λ * ΔF)

        # Logging
        msg = "Residual error is $err at iteration $iter."
        println(msg)
        open(error_file, "a") do file
            write(file, msg)
        end

        # Update history *after* using old values
        copy!(X_prev, X_k_original)
        copy!(F_prev, F_k)
        iter += 1
    end

    if err > tole
        error("PicardResidual failed to converge after $max_iters iterations. Final error = $err")
    end

    println()
    return solution_variables
end