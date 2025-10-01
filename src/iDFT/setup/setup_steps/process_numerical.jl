"""
    process_numerical(dataset::Dict) -> Dict{String, Float64}

Extracts and processes the `"numerical"` section from the dataset if it exists.
Returns a dictionary of numerical parameters as `Dict{String, Float64}`.
If no `"numerical"` section is present, returns an empty dictionary.
"""
function process_numerical(dataset::Dict)
    if haskey(dataset, "numerical")
        raw = dataset["numerical"]
        return Dict(string(k) => Float64(v) for (k, v) in raw)
    else
        return Dict{String, Float64}()
    end
end
