# === External Field Dispatchers ===
# These define cartesian_features(features, ::Val{:field_type}) and
# update_features(features, ::Val{:field_type}) for each supported type.

include("hard_walls.jl")         # :hard_walls
include("charged_walls.jl")      # :charged_walls
include("pointcharge_walls.jl")  # :pointcharge_walls
include("sw_walls.jl")           # :sw_walls
include("input_external.jl")     # :input_external
