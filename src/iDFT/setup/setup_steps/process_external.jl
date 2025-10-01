include("../../coordinate_system/coord_sys.jl")

function process_external(molsys, geometry)
    external_field = get!(geometry.features, :external_field, Dict{Any, Any}())

    @. geometry.features[:mirrored] = true
    
    species_properties = molsys.properties.species

    max_value = length(species_properties[:monomers])
    max_dim = length(geometry.features[:dimensions])

    for (_, field) in external_field
        if haskey(field, "monomers")
            for (i, monomer) in enumerate(species_properties[:monomers])
                monomer_str = string(monomer)  # Convert symbol to string to match YAML keys

                if haskey(field["monomers"], monomer_str)
                    monomer_props = field["monomers"][monomer_str]

                    for (property_str, value) in monomer_props
                        property_sym = Symbol(property_str)
                        surface_key = Symbol("surface_", property_str, "s")

                        if !haskey(field, surface_key)
                            field[surface_key] = [ [ zeros(Float64, 2) for _ in 1:max_value ] for _ in 1:max_dim ]
                        end

                        if isa(value, Vector{Vector{Float64}})
                            for (k, dim) in enumerate(field["dims"])
                                field[surface_key][dim][i] = value[k]
                            end
                        else
                            error("Expected property '$property_str' to be Vector{Vector{Float64}}.")
                        end
                    end
                end
            end
        end
    end
end