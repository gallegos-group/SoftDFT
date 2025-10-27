# === External Field Feature Dispatch ===
# Dispatches per-feature setup logic during grid preprocessing

include("external_fields/external_fields.jl")

"""
    cartesian_features(features::Dict)

Calls the appropriate `cartesian_features(features, ::Val{:feature})` method for each
entry in `features[:external_field]`. This applies external field metadata such as
setting `periodic = false`, storing wall positions, charges, etc.
"""
function cartesian_features(features::CartesianFeatures, data_external_field)
    external_field = Vector{AbstractExternalField}()
    for (name, _) in data_external_field
        push!(external_field, cartesian_features(features, data_external_field, Val{Symbol(name)}()))
    end

    return external_field
end

"""
    cartesian_features(features, ::Val{F})

Fallback method called when a specific handler is not defined.
"""
function cartesian_features(features::CartesianFeatures, data_external_field, ::Val{F}) where F
    error("No feature handler defined for external field: $F")
end

"""
    update_features(features::Dict)

Calls the appropriate `update_features(features, ::Val{:feature})` method for each
entry in `features[:external_field]`. Typically used to apply mirrored-domain corrections
(e.g., double-counting surface charges).
"""
function update_features(features::CartesianFeatures, data_external_field)
    for (name, _) in data_external_field
        update_features(features, data_external_field, Val{Symbol(name)}())
    end
end

"""
    update_features(features, ::Val{F})

Fallback method called when a specific update handler is not defined.
"""
function update_features(features::CartesianFeatures, data_external_field, ::Val{F}) where F
    error("No update handler defined for external field: $F")
end



# # === External Field Feature Dispatch ===
# # Dispatches per-feature setup logic during grid preprocessing

# include("external_fields/external_fields.jl")

# """
#     cartesian_features(features::Dict)

# Calls the appropriate `cartesian_features(features, ::Val{:feature})` method for each
# entry in `features[:external_field]`. This applies external field metadata such as
# setting `periodic = false`, storing wall positions, charges, etc.
# """
# function cartesian_features(features::Dict{Symbol, Any})
#     for external_field in features.external_field
#         cartesian_features(features, external_field)
#     end
# end

# """
#     cartesian_features(features, ::Val{F})

# Fallback method called when a specific handler is not defined.
# """
# function cartesian_features(features::Dict, external_field) where F
#     error("No feature handler defined for external field: $external_field")
# end

# """
#     update_features(features::Dict)

# Calls the appropriate `update_features(features, ::Val{:feature})` method for each
# entry in `features[:external_field]`. Typically used to apply mirrored-domain corrections
# (e.g., double-counting surface charges).
# """
# function update_features(features::Dict{Symbol, Any})
#     for external_field in features.external_field
#         update_features(features, external_field)
#     end
# end

# """
#     update_features(features, ::Val{F})

# Fallback method called when a specific update handler is not defined.
# """
# function update_features(features::Dict, external_field) where F
#     error("No update handler defined for external field: $external_field")
# end
