# === External Field Feature Dispatch ===
# Dispatches per-feature setup logic during grid preprocessing

include("external_fields/external_fields.jl")

"""
    cartesian_features(features::Dict)

Calls the appropriate `cartesian_features(features, ::Val{:feature})` method for each
entry in `features[:external_field]`. This applies external field metadata such as
setting `periodic = false`, storing wall positions, charges, etc.
"""
function cartesian_features(features::Dict{Symbol, Any})
    for (name, _) in features[:external_field]
        cartesian_features(features, Val{Symbol(name)}())
    end
end

"""
    cartesian_features(features, ::Val{F})

Fallback method called when a specific handler is not defined.
"""
function cartesian_features(features::Dict, ::Val{F}) where F
    error("No feature handler defined for external field: $F")
end

"""
    update_features(features::Dict)

Calls the appropriate `update_features(features, ::Val{:feature})` method for each
entry in `features[:external_field]`. Typically used to apply mirrored-domain corrections
(e.g., double-counting surface charges).
"""
function update_features(features::Dict{Symbol, Any})
    for (name, _) in features[:external_field]
        update_features(features, Val{Symbol(name)}())
    end
end

"""
    update_features(features, ::Val{F})

Fallback method called when a specific update handler is not defined.
"""
function update_features(features::Dict, ::Val{F}) where F
    error("No update handler defined for external field: $F")
end
