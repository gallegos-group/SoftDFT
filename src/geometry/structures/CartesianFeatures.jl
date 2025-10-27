struct CartesianFeatures
    dimensions    :: Vector{Float64}
    periodic      :: Vector{Bool}
    mirrored      :: Vector{Bool}
    offset        :: Vector{Int}
    total_charge  :: Vector{Float64}
end
