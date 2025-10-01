"""
    process_monomers(used_monomers, dataset, properties_struct)

Populates monomer-level physical properties (e.g., diameters, valences) for each monomer in `used_monomers`
based on the `"monomers"` section of the input dataset.

### Arguments
- `used_monomers::Dict{String, Int}`: Mapping from monomer names to their global indices.
- `dataset::Dict`: Parsed YAML input containing the `"monomers"` field.
- `properties_struct::PropertiesStruct`: Target structure where the `monomers` field will be filled.

### Behavior
- Iterates over all monomers in `used_monomers`, verifying each exists in the dataset.
- For each monomer, assigns its physical properties to the correct index range in `properties_struct.monomers`.
- Scalar properties are stored under pluralized keys (e.g., `"diameters"`, `"valences"`).
- Nested properties (e.g., interaction strengths with other monomers) are stored under keys like `"attraction: Y*s"`.

### Example
A YAML input like:
```yaml
monomers:
  A:
    diameter: 1.0
    valence: 1
    attraction:
      Y*: 2.0
Produces entries such as:
monomers["diameters"]         = [...]
monomers["valences"]          = [...]
monomers["attraction: Y*s"]   = [...]
"""
function process_monomers(used_monomers::Dict{String, Int}, dataset::Dict, properties_struct::PropertiesStruct)
    monomer_data = dataset["monomers"]
    monomer_properties = properties_struct.monomers

    max_index = maximum(values(used_monomers))
    sorted_monomers = sort(collect(used_monomers), by = x -> x[2])  # [(species, index)]

    prev_index = 0
    for (species, index) in sorted_monomers
        if !haskey(monomer_data, species)
            error("Species $species listed in used_monomers not found in dataset[\"monomers\"]")
        end

        props = monomer_data[species]
        for (property_str, value) in props
            prop = Symbol(property_str * "s")

            if isa(value, Dict)
                # Handle nested dictionary property (e.g., pKa: {acid: 4.5, base: 10.2})
                if !haskey(monomer_properties, prop)
                    monomer_properties[prop] = Dict{Symbol, Vector{Float64}}()
                end
                for (subkey_str, subval) in value
                    subkey = Symbol(subkey_str * "s")
                    if !haskey(monomer_properties[prop], subkey)
                        monomer_properties[prop][subkey] = zeros(max_index)
                    end
                    monomer_properties[prop][subkey][prev_index+1:index] .= subval
                end
            else
                # Handle flat scalar property (e.g., diameter: 1.0)
                if !haskey(monomer_properties, prop)
                    monomer_properties[prop] = zeros(max_index)
                end
                monomer_properties[prop][prev_index+1:index] .= value
            end
        end

        prev_index = index
    end
end


# function process_monomers(used_monomers::Dict{String, Int}, dataset::Dict, properties_struct::PropertiesStruct)
#     monomer_data = dataset["monomers"]
#     monomer_properties = properties_struct.monomers

#     max_index = maximum(values(used_monomers))
#     sorted_monomers = sort(collect(used_monomers), by = x -> x[2])  # [(species, index)]

#     prev_index = 0
#     for (species, index) in sorted_monomers
#         if !haskey(monomer_data, species)
#             error("Species $species listed in used_monomers not found in dataset[\"monomers\"]")
#         end

#         props = monomer_data[species]
#         for (property, value) in props
#             propertys = property * "s"  # plural key like "diameters", "valences"

#             if isa(value, Dict)
#                 # Property is a dict of subproperties
#                 for (subkey, subval) in value
#                     subkey_str = "$property: $subkey" * "s"
#                     if !haskey(monomer_properties, subkey_str)
#                         monomer_properties[subkey_str] = zeros(Float64, max_index)
#                     end
#                     monomer_properties[subkey_str][prev_index+1:index] .= subval
#                 end
#             else
#                 # Property is a scalar (e.g. diameter, valence)
#                 if !haskey(monomer_properties, propertys)
#                     monomer_properties[propertys] = zeros(Float64, max_index)
#                 end
#                 monomer_properties[propertys][prev_index+1:index] .= value
#             end
#         end

#         prev_index = index
#     end
# end
