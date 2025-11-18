# -------------------------------------------------------------------------
# functional_registry.jl
#
# This module defines the registry system for managing free energy 
# *functionals* used in inhomogeneous DFT (e.g., ideal gas, MFMT, MSA).
#
# Each functional contributes to the total Helmholtz free energy of the 
# system and is implemented as a callable object with an associated 
# constructor (typically `construct_functional(::Type{<:FreeEnergyFunctional}, ...)`).
#
# Key Functions:
# - `register_functional(name::Symbol, constructor::Function)`:
#       Registers a functional constructor under a given symbol key.
# - `get_functional(name::Symbol)`:
#       Retrieves a constructor by key, or throws an error if unknown.
# - `validate_functional_keys(keys::Vector{Symbol})`:
#       Checks that all provided keys exist in the registry.
# - `list_available_functionals()`:
#       Returns a sorted list of all registered functional keys.
#
# Usage example:
#     fn_type = get_functional(:mfmt)
#     fn_struct = construct_functional(fn_type, molsys, geometry)
# -------------------------------------------------------------------------

const FunctionalRegistry = Dict{Symbol, Type{<:AbstractFunctional}}()

function register_functional(name::Symbol, functional_type::Type{<:AbstractFunctional})
    FunctionalRegistry[name] = functional_type
end

function get_functional(name::Symbol)::Type{<:AbstractFunctional}
    get(FunctionalRegistry, name) do
        error("Unknown functional key: $name")
    end
end

function validate_functional_keys(keys::Vector{Symbol})
    unknown = filter(k -> !haskey(FunctionalRegistry, k), keys)
    if !isempty(unknown)
        error("Unknown functional(s): " * join(string.(unknown), ", "))
    end
end

function list_available_functionals()
    sort(collect(keys(FunctionalRegistry)))
end


# -------------------------------
# Explicit functional registrations
# -------------------------------

register_functional(:ideal, IdealFunctional)
# Hard-sphere repulsion
include("../excluded_volume/mfmt.jl")
register_functional(:mfmt, mfmtFunctional)

include("../excluded_volume/TPT1mfmt.jl")
register_functional(:TPT1mfmt, TPT1mfmtFunctional)

# Electrostatic correlations (MSA)
include("../electrostaticcorrelations/msa.jl")
register_functional(:msa, msaFunctional)

include("../electrostaticcorrelations/TPT1msa.jl")
register_functional(:TPT1msa, TPT1msaFunctional)

# Square-well interactions
include("../short_range_interactions/square_well.jl")
register_functional(:sw, SquareWellFunctional)