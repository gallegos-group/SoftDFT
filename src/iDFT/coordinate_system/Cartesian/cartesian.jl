# === Cartesian Coordinate System: Evaluation Tools ===

# External field evaluations for walls, charges, and user input
include("external_fields/external_fields.jl")

# Fourier transforms for common functions used in convolutions (step, delta, Coulomb)
include("process_fourier.jl")
include("fourier_transforms/fourier_transforms.jl")

# Fixed species
include("fixed_species/determine_fixed.jl")

# Fixed species
include("species_model/species_model.jl")

# Utilities for mirroring, domain halving, etc.
include("utilities/utilities.jl")

