function eval_rho_hat!(bulk_system :: IsingLST, geometry :: CartesianCoord, fields :: SpatialFields)
    
    calc_rho_beads_hat!(geometry, fields)

    calc_rho_bonds_hat!(bulk_system, geometry, fields)
end