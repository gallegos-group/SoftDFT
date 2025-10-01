# -------------------------------------------------------------------------
# AbstractFreeEnergy.jl
#
# This module defines abstract types for thermodynamic model components
# used in bulk theories (e.g., classical DFT, liquid-state theory, and
# chain-based theories).
#
# Each model contributes a specific Helmholtz free energy term such as 
# ideal entropy, excluded volume, Coulombic interactions, or bonding.
#
# Abstract types enable consistent typing, dispatch, and trait logic.
#
# Key Abstract Types:
# - `AbstractFreeEnergy`: Base class for all bulk thermodynamic models.
# - `IdealFreeEnergy`:    Subclass for ideal-gas contributions.
# - `ExcessFreeEnergy`:   Subclass for all non-ideal contributions.
#
# Trait Utilities:
# - `is_ideal(model::AbstractFreeEnergy)`: returns true for ideal terms.
#
# Note: Model registration is handled externally in `model_registry.jl`.
# -------------------------------------------------------------------------

# -------------------------
# Abstract Type Hierarchy
# -------------------------
abstract type AbstractFreeEnergy end

abstract type ExcessFreeEnergy <: AbstractFreeEnergy end

# -------------------------
# Trait Dispatch
# -------------------------
"""
    is_ideal(model::AbstractFreeEnergy) -> Bool

Trait function that returns `true` for ideal gas models,
`false` for all others by default.
"""
is_ideal(::AbstractFreeEnergy) = false
