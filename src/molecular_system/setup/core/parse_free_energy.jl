"""
    parse_(models_list_raw)

Convert the list of user-specified model names into standardized `Symbol` identifiers.

Given:
- A single model name (`String`) or a list of model names (`Vector{String}`),

This function:
1. Ensures the input is treated as a list of strings,
2. Converts each string into a `Symbol`,
3. Returns a list of model symbols suitable for internal use and validation.

Returns:
- `models_list::Vector{Symbol}` — list of model names as symbols, e.g., `[:ideal, :sw, :msa]`.

Example:
```julia
process_model("ideal")              # → [:ideal]
process_model(["ideal", "msa"])    # → [:ideal, :msa]
"""

function parse_free_energy(models_list_raw)
    model_strings = isa(models_list_raw, String) ? [models_list_raw] : models_list_raw

    models_list = Symbol[]
    for name in model_strings
        push!(models_list, Symbol(name))
    end

    return models_list
end