"""
    get_temperature(dataset::Dict) -> Float64

Extract the temperature (in Kelvin) from the input dataset.

If the value is specified under `constraints["temperature"]`, it is returned after validation.
If not provided, defaults to 298.15 K (room temperature).

Throws an error if the specified value is not a `Float64`.
"""
function get_temperature(dataset)
    temperature = 298.15

    if haskey(dataset, "constraints")
        if haskey(dataset["constraints"], "temperature")
            temperature_value = dataset["constraints"]["temperature"]
            temperature = temperature_value isa Float64 ? temperature_value : error("Temperature must be a Float64 value.")
        else
            println("No temperature specified. Defaulting to 298.15 K")
        end
    end

    return temperature
end