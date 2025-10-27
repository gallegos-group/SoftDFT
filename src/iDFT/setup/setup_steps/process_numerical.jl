"""
    process_numerical(dataset::Dict) -> Dict{String, Float64}

Extracts and processes the `"numerical"` section from the dataset if it exists.
Returns a dictionary of numerical parameters as `Dict{String, Float64}`.
If no `"numerical"` section is present, returns an empty dictionary.
"""
# function process_numerical(dataset::Dict)
#     if haskey(dataset, "numerical")
#         raw = dataset["numerical"]
#         return Dict(string(k) => Float64(v) for (k, v) in raw)
#     else
#         return Dict{String, Float64}()
#     end
# end

function process_numerical(dataset::Dict{Any, Any})
    if haskey(dataset, "numerical")
        raw = dataset["numerical"]
        out = Dict{String, Float64}()
        @inbounds for (k, v) in raw
            out[k] = v isa Real ? Float64(v) : parse(Float64, string(v))
        end
        return out
    else
        return Dict{String, Float64}()
    end
end

