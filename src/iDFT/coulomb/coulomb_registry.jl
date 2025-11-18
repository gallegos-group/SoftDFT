# Coulomb Model Registry
# Maps symbol → Coulomb model type

const coulomb_registry = Dict{Symbol, Type}(
    :coul => CoulombFFT,
)

"""
    get_coulomb_model(name)

Return the Coulomb model type corresponding to a symbol or string `name`.
"""
function get_coulomb_model(name)
    key = name isa Symbol ? name : Symbol(name)

    if haskey(coulomb_registry, key)
        return coulomb_registry[key]
    end
end
