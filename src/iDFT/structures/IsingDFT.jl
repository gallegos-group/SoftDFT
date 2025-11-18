include("SpatialFields.jl")
include("../coulomb/AbstractCoul.jl")
include("../functionals/Abstract_Types/AbstractFunctional.jl")

"""
    IsingDFT(bulk_system, geometry, fields, functionals, coulomb, numerics)

Top-level container for an inhomogeneous Ising DFT calculation.

This struct ties together:
- the molecular model and bulk thermodynamics (`bulk_system`),
- the spatial discretization (`geometry`),
- all iteration-dependent fields (`fields`),
- local excess free-energy functionals (`functionals`),
- the electrostatic model (`coulomb`),
- and numerical solver controls (`numerics`).

# Type Parameters
- `TIsingLST` : Type of the bulk Ising liquid-state theory object
  (e.g. `IsingLST`), which holds the molecular system (`molsys`) and
  bulk thermodynamic state.
- `TGeom`     : Coordinate/geometry type, typically `CartesianGrid{D}`.
- `TFields`   : Spatial fields container, typically `SpatialFields`,
  holding densities, excess fields (`mu_ex_K`, `lng_K`, `Psi`, `PsiC`),
  Fourier caches, and fixed species data.
- `TFuncs`    : Tuple of **local** excess free-energy functionals
  (e.g. mFMT, TPT1, association, short-range interactions). Coulomb
  is intentionally excluded from this tuple.
- `TCoul`     : Coulomb / electrostatic model (e.g. a `coulFunctional`
  containing the k-space kernel and electrostatic parameters).

# Fields
- `bulk_system :: TIsingLST`
    Molecular system and bulk thermodynamic problem. Includes species,
    sequences, state families, and bulk reference state.
- `geometry    :: TGeom`
    Spatial discretization and boundary conditions (grid dimensions,
    bin widths, periodic/mirrored flags, offsets, total fixed charge,
    and geometric features).
- `fields      :: TFields`
    Iteration-dependent spatial fields:
      * densities (`rho_K`),
      * excess fields (`lambda_K`, `mu_ex_K`, `lng_K`, `Ext`, `Psi`, `PsiC`),
      * Fourier-space work arrays,
      * fixed/grafted segment information.
- `functionals :: TFuncs`
    Tuple of local (short-range) excess free-energy functionals that
    contribute to `mu_ex_K` and `lng_K`. Evaluated via
    `eval_mu_functionals` and **do not** include long-range Coulomb.
- `coulomb     :: TCoul`
    Electrostatic model used to build and solve the Poisson equation
    (e.g. via FFT using precomputed `weight_fft_K`) and to determine
    the neutrality correction (`PsiC`).
- `numerics    :: Dict{String, Float64}`
    Numerical parameters for the nonlinear solver (tolerances, maximum
    iterations, mixing parameters, step sizes, etc.).

# Usage
Construct an `IsingDFT` instance once, then pass it to a solver
(e.g. `new_solver_dft`) which will:
1. Update densities using `fields.excess.mu_ex_K`, `Psi`, and `PsiC`,
2. Evaluate local excess functionals (`functionals`),
3. Solve the Poisson equation using `coulomb` and `geometry`,
4. Iterate on `(mu_ex_K, lng_K, Psi, PsiC)` using the settings in
   `numerics` until self-consistency is reached.
"""
struct IsingDFT{TIsingLST, TGeom, TFields, TFuncs, TCoul}
    bulk_system :: TIsingLST
    geometry    :: TGeom
    fields      :: TFields
    functionals :: TFuncs
    coulomb     :: TCoul
    numerics    :: Dict{String, Float64}
end
