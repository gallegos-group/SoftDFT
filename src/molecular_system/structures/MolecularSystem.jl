"""
    MolecularSystem

Top-level data structure representing a fully specified molecular system,
used in bulk or inhomogeneous simulations. It encapsulates all polymer configurations,
chemical species, reaction definitions, model selections, and physical constants.

# Fields
- `configurations::Vector{ConfigurationStruct}`: List of polymer or molecular species
  to simulate. Each has its own sequence, topology, states, and associated physics model.

- `properties::PropertiesStruct`: Contains parsed chemical definitions, reaction networks,
  monomer-level properties, and selected free energy models. Used during density and
  energy evaluations.

- `constants::ConstantsStruct`: Physical constants (e.g., `k_B`, `e₀`, `N_A`) and
  user-defined length scale (`l_sc`). Used to compute derived quantities such as
  the Bjerrum length.

This struct is constructed by `setup_molsys()` and passed to all downstream components,
including solvers, density functionals, and thermodynamic analysis routines. It serves as
the canonical system definition.
"""


# These define the foundational types used throughout the molecular system.
include("ConstantsStruct.jl")       # Physical constants (k_B, e0, etc.)
include("ReactionStruct.jl")        # Chemical reaction definitions
include("TopologyStruct.jl")        # Topology
include("ConfigurationStruct.jl")   # Per-species polymer model and sequence info
include("PropertiesStruct.jl")      # Species/reaction properties and FE models

struct MolecularSystem{CS}
    configurations::CS
    properties::PropertiesStruct
    constants::ConstantsStruct
end
