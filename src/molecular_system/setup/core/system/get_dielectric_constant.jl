"""
    get_dielectric_constant(dataset::Dict) -> Float64

Extract the relative dielectric constant (ε_r) from the input dataset.

If the value is specified under `constraints["dielectric_constant"]`, it is returned after validation.
If not provided, defaults to 78.4 (corresponding to water at 25°C).

Throws an error if the specified value is not a `Float64`.
"""
function get_dielectric_constant(dataset)
    
    dielectric_constant = 78.4

    if haskey(dataset, "constraints")
        if haskey(dataset["constraints"], "dielectric_constant")
            dielectric_value = dataset["constraints"]["dielectric_constant"]
            dielectric_constant = dielectric_value isa Float64 ? dielectric_value : error("Dielectric constant must be a Float64 value.")
        else
            println("No dielectric value specified. Defaulting to 78.4.")
        end
    end

    return dielectric_constant
end
