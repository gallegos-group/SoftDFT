"""
    freely_jointed.jl

Includes all function definitions related to the `freely_jointed` species model 
used in bulk density evaluations for chain-like molecules.

This model assumes a freely-jointed chain architecture, enabling recursive 
propagator calculations (forward and backward) over a tree-like topology.

Components included:
- `child_propagators.jl`: Computes forward (child) propagators `gC`.
- `parent_propagators.jl`: Computes backward (parent) propagators `gP`.
- `partition_function.jl`: Evaluates the partition function `Xi` for normalization.
- `accumulate_densities.jl`: Computes bead and bond density contributions using propagators.
- `eval_bulk_density.jl`: Orchestrates propagator setup, partition function evaluation,
  and density accumulation for this species model.

This file is meant to be included once in the top-level solver/module setup:
    include("species_model/freely_jointed.jl")
"""

include("freely_jointed/child_propagators.jl")
include("freely_jointed/parent_propagators.jl")
include("freely_jointed/partition_function.jl")
include("freely_jointed/accumulate_densities.jl")
include("freely_jointed/eval_bulk_density.jl")
include("freely_jointed/bulk_chemical_potential.jl")