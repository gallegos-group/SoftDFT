"""
    calc_rho_beads_hat!(geometry, fields)

Compute Fourier–space bead densities from real–space bead densities.

For each trailing index `j`:
1. Fill `f_hat` with `rho_beads_K .* trapez` on the true domain.
2. Map into storage indices and apply mirroring.
3. Forward FFT in-place (`f_hat = plan_forward * f_hat`).
4. Store result in `rho_beads_hat`.

Both `f_hat` and `rho_beads_hat` are defined on the full
calculation domain `NP_star`.
"""


function calc_rho_beads_hat!(geometry :: CartesianCoord, fields :: SpatialFields)
    @unpack NP, mirrored, offset = geometry

    NP_star = compute_full_domain(NP, mirrored, offset)

    rho_beads_K = fields.rho_K.beads
    @unpack Ext, trapez = fields.excess
    @unpack rho_beads_hat, f_hat, plan_forward = fields.fourier

    Rsys = CartesianIndices(NP)
    Rsys_star = CartesianIndices(NP_star)

    for j in axes(rho_beads_K, ndims(rho_beads_K))
        @. f_hat = 0.0
        
        @inbounds for K in Rsys
            K_star = to_star_index(K, offset, mirrored)
            
            f_hat[K_star] = rho_beads_K[K, j] * trapez[K, j]
        end

        # extend to mirrored halves
        apply_mirroring!(f_hat, geometry)        

        # FFT in-place
        f_hat = plan_forward * f_hat

        @inbounds for K_star in Rsys_star
            rho_beads_hat[K_star, j] = f_hat[K_star]
        end
    end
end

"""
    calc_rho_bonds_hat!(bulk_system, geometry, fields)

Compute Fourier–space bond densities from real–space bond densities.

For each bond type `(i, j)` and bond index `u`, both bond ends
are transformed:
1. Fill `f_hat` with `rho_bonds_K .* trapez` for the selected end.
2. Apply mirroring into the full domain.
3. Forward FFT in-place (`f_hat = plan_forward * f_hat`).
4. Store result in `rho_bonds_hat`.

Both `f_hat` and `rho_bonds_hat` are defined on the full
calculation domain `NP_star`.
"""

function calc_rho_bonds_hat!(bulk_system :: IsingLST, geometry :: CartesianCoord, fields :: SpatialFields)

    @unpack NP, mirrored, offset = geometry

    NP_star = compute_full_domain(NP, mirrored, offset)

    rho_bonds_K = fields.rho_K.bonds
    trapez = fields.excess.trapez
    @unpack rho_bonds_hat, f_hat, plan_forward = fields.fourier

    Rsys      = CartesianIndices(NP)
    Rsys_star = CartesianIndices(NP_star)

    for (u, (i, j)) in enumerate(bulk_system.molsys.properties.bond_types)
        for (idx, bead_type) in ((1, i), (2, j))
            @. f_hat = 0.0

            @inbounds for K in Rsys
                K_star = to_star_index(K, offset, mirrored)
                f_hat[K_star] = rho_bonds_K[K, idx, u] * trapez[K, bead_type]
            end

            apply_mirroring!(f_hat, geometry)

            f_hat = plan_forward * f_hat

            @inbounds for K_star in Rsys_star
                rho_bonds_hat[K_star, idx, u] = f_hat[K_star]
            end
        end
    end
end