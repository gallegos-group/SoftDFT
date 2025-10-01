"""
    compute_bjerrum_length(constants::ConstantsStruct, εr::Float64, T::Float64) -> Float64

Compute the Bjerrum length in simulation units (scaled by `l_sc`):

    ℓ_B = e² / (4π ε₀ ε_r k_B T)

This is the distance at which two unit charges interact with thermal energy ( k_B T ).
The result is returned in units of `l_sc`.
"""

function compute_bjerrum_length(constants::ConstantsStruct, εr::Float64, T::Float64)
    return constants.e0^2 / (4 * π * constants.epsilon0 * εr * constants.k_B * T) / constants.l_sc
end
