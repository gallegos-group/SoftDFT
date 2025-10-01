"""
    process_system(dataset::Dict, constants::ConstantsStruct) -> Dict{String, Float64}

Parse and compute system-level physical properties from the input dataset.

Returns a dictionary with:
- `"temperature"`: System temperature in Kelvin
- `"dielectric_constant"`: Relative permittivity ε_r
- `"bjerrum_length"`: Bjerrum length in simulation units (scaled by `l_sc`)

The Bjerrum length is computed as:
    ℓ_B = e² / (4π ε₀ ε_r k_B T) / l_sc
"""
function process_system(dataset, constants::ConstantsStruct)
    system_properties = Dict{String, Float64}()

    T = get_temperature(dataset)
    εr = get_dielectric_constant(dataset)
    BJ = compute_bjerrum_length(constants, εr, T)

    system_properties["temperature"] = T
    system_properties["dielectric_constant"] = εr
    system_properties["bjerrum_length"] = BJ

    return system_properties
end
