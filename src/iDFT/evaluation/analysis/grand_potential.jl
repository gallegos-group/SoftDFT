include("excess_free_energy.jl")
include("num_species.jl")

function grand_potential(dft_system :: IsingDFT)

    @unpack bulk_system, geometry, fields, functionals = dft_system

    rho_beads_K = fields.rho_K.beads
    rho_bonds_K = fields.rho_K.bonds

    @unpack mu_ex_K, trapez, lng_K = fields.excess

    Omega = eval_fe_functional(bulk_system, geometry, fields, functionals)
    
    for Kj in CartesianIndices(rho_beads_K)
        Omega -= rho_beads_K[Kj] * mu_ex_K[Kj] * trapez[Kj]
    end
    
    for Kj in CartesianIndices(rho_bonds_K)
        Omega += rho_bonds_K[Kj] * lng_K[Kj] / 2.0
    end
 
    Omega *= prod(geometry.bin_width)

    Omega -= sum(get_num_species(bulk_system, geometry, fields))

    return Omega
end