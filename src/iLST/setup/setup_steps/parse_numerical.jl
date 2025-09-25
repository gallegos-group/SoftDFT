"""
    parse_numerical(numerical_data)

Converts a dictionary of numerical parameters (possibly with mixed types) into
a `Dict{String, Float64}`. All values are converted using `Float64(x)`.

Throws an error if any value cannot be safely converted.
"""

function parse_numerical(dataset)
    
    numerical_data = get!(dataset, "numerical", Dict{String, Float64}())
    numerics = Dict{String, Float64}()

    for (k, v) in numerical_data
        try
            numerics[k] = Float64(v)
        catch
            error("Could not convert value for key '$k' to Float64: $v")
        end
    end

    return numerics
end
