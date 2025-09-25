
function output_arrays(dft_system :: IsingDFT)
    write_field_array(dft_system.fields.excess.Psi, "MEP.txt", dft_system.geometry.NP; names=["MEP"])
    write_field_array(dft_system.fields.rho_K.beads, "rho_beads_K.txt", dft_system.geometry.NP; names=dft_system.bulk_system.molsys.properties.species[:monomers])
    write_segment_array(dft_system.fields.rho_K.segments, dft_system.bulk_system.molsys, dft_system.geometry.NP)
    write_state_fraction_profiles(dft_system.fields.rho_K.segments,dft_system.bulk_system.molsys,dft_system.geometry.NP)

    write_total_state_fractions(dft_system.fields.rho_K.segments,dft_system.bulk_system.molsys,dft_system.geometry.NP)
end