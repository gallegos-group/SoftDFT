# -------------------------------------------------------------------------
# AbstractFunctional.jl
#
# This module defines abstract types for representing Helmholtz free energy
# *functionals* used in inhomogeneous DFT (iDFT) simulations.
#
# Each functional contributes a specific term to the total free energy 
# (e.g., ideal gas, hard-sphere repulsion, electrostatics, bonding, etc.)
# and must be implemented as a subtype of `AbstractFunctional`.
#
# Key Abstract Types:
# - `AbstractFunctional`: Base class for all iDFT free energy functionals.
# - `IdealFunctional`:    For ideal gas contributions.
# - `ExcessFunctional`:   For all non-ideal terms (e.g. MFMT, Coulomb, MSA).
#
# Optional traits (e.g., `is_ideal`) can be used for specialized logic.
#
# Usage:
#     struct MyFunctional <: ExcessFunctional ... end
#     construct_functional(::Type{MyFunctional}, molsys, geom) = ...
# -------------------------------------------------------------------------

# -------------------------
# Abstract Type Hierarchy
# -------------------------
abstract type AbstractFunctional end

abstract type IdealFunctional  <: AbstractFunctional end
abstract type ExcessFunctional <: AbstractFunctional end

# -------------------------
# Functional Trait Utilities
# -------------------------
"""
    is_ideal(f::AbstractFunctional) -> Bool

Returns true if the functional is an ideal gas contribution.
Default is `false`; overridden for `IdealFunctional` types.
"""
is_ideal(::Type{<:AbstractFunctional}) = false
is_ideal(::Type{IdealFunctional})    = true
