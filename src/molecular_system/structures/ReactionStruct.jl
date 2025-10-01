"""
    ReactionStruct

Defines a chemical reaction between one or more species, typically used to model
protonation, deprotonation, or other association/dissociation equilibria.
Each reaction is defined by its species participants, stoichiometric coefficients,
and equilibrium constant.

Fields:
- `species::Vector{String}`: Names of the species involved in the reaction, e.g. `["A", "H", "B"]`.
  Species can be monomers, ions, or neutral molecules.

- `coeffs::Vector{Int}`: Stoichiometric coefficients corresponding to each species.
  Negative values represent reactants, positive values represent products.
  Must be ±1 to ensure well-defined logarithmic mass-action expressions.

- `pK::Vector{Float64}`: Base-10 logarithm of the equilibrium constant for the reaction,
  i.e. `pK = -log10(K_eq)`. Only one value is expected per reaction.

This struct is used during molecular system setup to construct state families,
evaluate chemical potentials, and enforce mass-action constraints during
self-consistent field iterations.
"""

struct ReactionStruct
  species  :: Vector{String}       # e.g., ["A", "H", "B"]
  coeffs   :: Vector{Int}          # e.g., [-1, -1, +1]
  pK       :: Vector{Float64}      # equilibrium constant in -log scale
end