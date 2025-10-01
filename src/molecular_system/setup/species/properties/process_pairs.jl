"""
    process_pairs(used_monomers, dataset, properties_struct)

Parses pairwise interaction properties between monomers (e.g., square-well energy, lambda values) 
from the `"pairs"` section of the input dataset and populates the `pairs` field in `properties_struct`.

### Arguments
- `used_monomers::Dict{String, Int}`: Mapping of monomer names to global indices.
- `dataset::Dict`: Parsed YAML input.
- `properties_struct::PropertiesStruct`: Target structure to populate with pair interaction data.

### Behavior
- For each defined pair type (e.g., `"square_well"`), generates:
  - `"square_well pairs"`: a list of unique `(i,j)` index pairs.
  - `"square_well energys"`: corresponding interaction energies.
  - `"square_well lambdas"`: interaction range multipliers (defaulted to `1.5` if not specified).
- Only pairs explicitly defined in the YAML are considered.
- Also populates `"square_well species"` with the list of monomer indices involved in any defined pair.

### Example
A YAML input like:
```yaml
pairs:
  square_well:
    "A, B":
      energy: 2.0
      lambda: 1.5
Creates:
pairs["square_well pairs"]    = [(1, 2)]
pairs["square_well energys"]  = [2.0]
pairs["square_well lambdas"]  = [1.5]
pairs["square_well species"]  = [1, 2]
"""

function process_pairs(used_monomers::Dict{String, Int}, dataset::Dict, properties_struct)
    pair_data = get(dataset, "pairs", Dict{String, Any}())
    pair_properties = properties_struct.pairs

    sorted_monomers = sort(collect(used_monomers), by = x -> x[2])  # [(A, 1), (B, 2), ...]

    for (pair_type_str, pairs_dict) in pair_data
        pair_type = Symbol(pair_type_str)  # e.g., :square_well, :sticky
        pair_properties[pair_type] = Dict{Symbol, Any}()

        pair_properties[pair_type][:pairs] = Tuple{Int, Int}[]

        for ((name1, index1), (name2, index2)) in Iterators.product(sorted_monomers, sorted_monomers)
            pair_name = string(name1, ", ", name2)

            if haskey(pairs_dict, pair_name)
                interaction = pairs_dict[pair_name]

                j1, j2 = index1, index2
                ordered_pair = j1 < j2 ? (j1, j2) : (j2, j1)
                push!(pair_properties[pair_type][:pairs], ordered_pair)

                for (property_str, value) in interaction
                    prop_key = Symbol(property_str * "s")  # standardize to plural field
                    if !haskey(pair_properties[pair_type], prop_key)
                        pair_properties[pair_type][prop_key] = Float64[]
                    end
                    push!(pair_properties[pair_type][prop_key], value)
                end
            end
        end

        # Post-processing: add :lambdas default if :energys exists but no :lambdas
        if :energys in keys(pair_properties[pair_type])
            n = length(pair_properties[pair_type][:energys])
            if !haskey(pair_properties[pair_type], :lambdas)
                pair_properties[pair_type][:lambdas] = fill(1.5, n)
            end
        end

        # Record which species are involved
        involved = Set{Int}()
        for (j1, j2) in pair_properties[pair_type][:pairs]
            push!(involved, j1)
            push!(involved, j2)
        end
        pair_properties[pair_type][:species] = collect(involved)
    end

    return pair_properties
end


# function process_pairs(used_monomers::Dict{String, Int}, dataset::Dict, properties_struct)
#     pair_data = get(dataset, "pairs", Dict{String, Any}())
#     pair_properties = properties_struct.pairs

#     sorted_monomers = sort(collect(used_monomers), by = x -> x[2])  # [(A, 1), (B, 2), ...]

#     for (pair_type, pairs_dict) in pair_data
#         pair_key = pair_type * " pairs"
#         property_prefix = pair_type * " "

#         pair_properties[pair_key] = Tuple{Int, Int}[]

#         for ((name1, index1), (name2, index2)) in Iterators.product(sorted_monomers, sorted_monomers)
#             pair_name = string(name1, ", ", name2)

#             if haskey(pairs_dict, pair_name)
#                 interaction = pairs_dict[pair_name]

#                 j1, j2 = index1, index2
#                 ordered_pair = j1 < j2 ? (j1, j2) : (j2, j1)
#                 push!(pair_properties[pair_key], ordered_pair)

#                 for (property, value) in interaction
#                     prop_name = property_prefix * property * "s"
#                     if !haskey(pair_properties, prop_name)
#                         pair_properties[prop_name] = Float64[]
#                     end
#                     push!(pair_properties[prop_name], value)
#                 end
#             end
#         end
#     end

#     # === Postprocessing Defaults ===
#     if haskey(pair_properties, "square_well energys")
#         # initialize lambdas if missing
#         n_pairs = length(pair_properties["square_well energys"])
#         if !haskey(pair_properties, "square_well lambdas")
#             pair_properties["square_well lambdas"] = fill(1.5, n_pairs)
#             # println("Note: Initialized missing 'square_well lambdas' to 1.5")
#         end

#         # track all unique species participating
#         involved = Set{Int}()
#         for (j1, j2) in pair_properties["square_well pairs"]
#             push!(involved, j1)
#             push!(involved, j2)
#         end
#         pair_properties["square_well species"] = collect(involved)
#     end

#     return pair_properties
# end
