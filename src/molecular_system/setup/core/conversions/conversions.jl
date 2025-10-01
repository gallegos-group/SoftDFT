# useful conversions

"""
    density_to_concentration(rho::Float64, constants::ConstantsStruct) -> Float64

Convert number density in units of 1/l_sc³ to mol/L.
"""
function density_to_concentration(rho::Float64, constants::ConstantsStruct)
    # Volume of 1 l_sc³ in liters: (l_sc [m])³ × 1000 L/m³
    vol_liters = (constants.l_sc)^3 * 1000.0  # m³ × 1000 = L
    return rho / vol_liters / constants.NA  # mol/L
end

"""
    concentration_to_density(conc::Float64, constants::ConstantsStruct) -> Float64

Convert concentration in mol/L to number density in units of 1/l_sc³.
"""
function concentration_to_density(conc::Float64, constants::ConstantsStruct)
    vol_liters = (constants.l_sc)^3 * 1000.0
    return conc * constants.NA * vol_liters
end


"""
    kBT_to_J(T::Float64, constants::ConstantsStruct) -> Float64

Returns k_BT in joules for a given temperature.
"""
kBT_to_J(T::Float64, constants::ConstantsStruct) = constants.k_B * T

"""
    J_to_kBT(E::Float64, T::Float64, constants::ConstantsStruct) -> Float64

Converts energy from joules to k_BT.
"""
J_to_kBT(E::Float64, T::Float64, constants::ConstantsStruct) = E / (constants.k_B * T)


"""
    mM_to_molL(c_mM) -> Float64

Convert concentration from millimolar (mM) to mol/L.
"""
mM_to_molL(c_mM::Float64) = c_mM * 1e-3

"""
    molL_to_mM(c_molL) -> Float64

Convert concentration from mol/L to millimolar (mM).
"""
molL_to_mM(c_molL::Float64) = c_molL * 1e3

"""
    kcalmol_to_J(E_kcalmol) -> Float64

Convert energy from kcal/mol to joules.
"""
kcalmol_to_J(E_kcalmol::Float64) = E_kcalmol * 4184.0 / 6.02214076e23  # exact NA

"""
    J_to_kcalmol(E_J) -> Float64

Convert energy from joules to kcal/mol.
"""
J_to_kcalmol(E_J::Float64) = E_J * 6.02214076e23 / 4184.0

"""
    nm_to_lsc(length_nm, constants) -> Float64

Convert a length from nanometers to reduced units (multiples of `l_sc`).
"""
nm_to_lsc(length_nm::Float64, constants::ConstantsStruct) = (length_nm * 1e-9) / constants.l_sc

"""
    lsc_to_nm(length_lsc, constants) -> Float64

Convert a reduced length (in units of `l_sc`) to nanometers.
"""
lsc_to_nm(length_lsc::Float64, constants::ConstantsStruct) = length_lsc * constants.l_sc * 1e9
