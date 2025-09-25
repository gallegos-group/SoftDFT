"""
    monomerbead.jl

Includes all function definitions related to the `monomerbead` species model
used in bulk density evaluations for single-segment (monomeric) molecules
with discrete internal states.

This model treats each species as a non-interacting monomer, omitting any
bonded connectivity. Only local chemical potentials and internal states are 
considered in the density evaluation.

Components included:
- `partition_function.jl`: Evaluates the partition function `Xi` for normalization.
- `accumulate_bulk_densities.jl`: Computes bead density contributions using `Xi`.
- `eval_bulk_density.jl`: Orchestrates evaluation of chemical potentials, normalization,
  and accumulation for this species model.

This file is meant to be included once in the top-level solver/module setup:
    include("species_model/monomerbead.jl")
"""

include("monomerbead/partition_function.jl")
include("monomerbead/accumulate_densities.jl")
include("monomerbead/eval_bulk_density.jl")
include("monomerbead/bulk_chemical_potential.jl")