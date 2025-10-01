# core/SolverUtils.jl

# Shared helpers
include("helpers/count_total_entries.jl")
include("helpers/apply_mixing.jl")
include("helpers/flat_solution.jl")
include("helpers/apply_flat_solution.jl")
include("helpers/pursuing_target.jl")

"""
    norm_converged(r, tol) → Bool

Check if the norm of residual or update is below a tolerance.
"""
norm_converged(r::AbstractVector, tol::Real) = norm(r) < tol
norm_converged(r::Real, tol::Real) = abs(r) < tol

"""
    has_diverged(x) → Bool

Detect NaNs or Infs in solution.
"""
has_diverged(x) = any(isnan, x) || any(isinf, x)

"""
    print_iteration(k, norm_r)

Print basic iteration diagnostics.
"""
function print_iteration(k::Int, norm_r::Real)
    @info "Iter $k: residual norm = $norm_r"
end

"""
    maybe_stop(k, max_iter) → Bool

Returns true if max iteration is reached.
"""
maybe_stop(k::Int, max_iter::Int) = k >= max_iter
