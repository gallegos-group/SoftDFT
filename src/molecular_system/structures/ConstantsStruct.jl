"""
    ConstantsStruct

Holds fundamental physical constants used throughout the molecular system setup,
thermodynamic calculations, and electrostatic models. Values are in SI units.

Fields:
- `epsilon0::Float64`: Vacuum permittivity (F·m⁻¹)
- `e0::Float64`: Elementary charge (C)
- `k_B::Float64`: Boltzmann constant (J·K⁻¹)
- `NA::Float64`: Avogadro's number (mol⁻¹)
- `l_sc::Float64`: Default molecular length scale, e.g. monomer size (m)

These values are used to compute derived quantities like the Bjerrum length,
dimensionless energies, and unit conversions. Defaults are provided but can be
overridden as needed.
"""
Base.@kwdef struct ConstantsStruct
    epsilon0::Float64 = 8.854187817e-12
    e0::Float64 = 1.60217653e-19
    k_B::Float64 = 1.3806505e-23
    NA::Float64 = 6.0221415e23
    l_sc::Float64 = 0.5e-9
end
