# -------------------------------------------------------------------------
# model_registry.jl
#
# This module defines the abstract type `AbstractFreeEnergy` and provides
# a registry system for managing thermodynamic model components used in
# SoftDFT and related solvers (e.g., classical DFT, liquid-state theory,
# or chain-based theories).
#
# Each model corresponds to a specific contribution to the total Helmholtz
# free energy (e.g., ideal, excluded volume, electrostatic), and is 
# registered using a unique symbol key. This modular registry enables 
# dynamic model composition at runtime.
#
# Key Functions:
# - `register_model(name::Symbol, model_type::Type{<:AbstractFreeEnergy})`:
#       Registers a model type under a given name.
# - `get_model(name::Symbol)`:
#       Retrieves the registered model type by name.
# - `validate_model_keys(keys::Vector{Symbol})`:
#       Verifies that all provided keys exist in the registry.
# - `list_available_models()`:
#       Returns a sorted list of registered model keys.
#
# Example:
#     model_type = get_model(:ideal)
#     model_instance = construct_free_energy(model_type, rho)
# -------------------------------------------------------------------------

const ModelRegistry = Dict{Symbol, Type{<:AbstractFreeEnergy}}()

function register_model(name::Symbol, model_type::Type{<:AbstractFreeEnergy})
    ModelRegistry[name] = model_type
end

function get_model(name::Symbol)
    get(ModelRegistry, name) do
        error("Unknown model key: $name")
    end
end

function validate_model_keys(keys::Vector{Symbol})
    unknown = filter(k -> !haskey(ModelRegistry, k), keys)
    if !isempty(unknown)
        error("Unknown model(s): " * join(string.(unknown), ", "))
    end
end

function list_available_models()
    sort(collect(keys(ModelRegistry)))
end

# -------------------------------
# Explicit model registrations
# -------------------------------

# Ideal contribution
include("../ideal/ideal.jl")
register_model(:ideal, IdealFreeEnergy)

# Excluded volume contributions
include("../excluded_volume/mfmt.jl")
register_model(:mfmt, MFMTFreeEnergy)

include("../excluded_volume/TPT1mfmt.jl")
register_model(:TPT1mfmt, TPT1MFMTFreeEnergy)

# Electrostatic correlations
include("../electrostaticcorrelations/msa.jl")
register_model(:msa, MSAFreeEnergy)

include("../electrostaticcorrelations/TPT1msa.jl")
register_model(:TPT1msa, TPT1MSAFreeEnergy)

# Short-range interactions
include("../short_range_interactions/square_well.jl")
register_model(:sw, SquareWellFreeEnergy)
