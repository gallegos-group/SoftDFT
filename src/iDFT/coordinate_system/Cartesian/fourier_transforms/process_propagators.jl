"""
    process_propagator(bulk_system, geometry::CartesianCoord) -> Array{Float64, N}

Compute the FFT-ready propagator weights for each bond type in the system.

# Arguments
- `bulk_system` :: Object containing bond definitions and monomer diameters.
- `geometry`    :: Cartesian coordinate system (used to define bin width and grid shape).

# Returns
- `weight_bond_hat` :: Real-space bond weight field for each bond type, ready for FFT.
"""
function process_propagator(molsys::MolecularSystem, geometry::CartesianCoord)
    @unpack bin_width, NP, mirrored, offset = geometry

    NP_star = compute_full_domain(NP, mirrored, offset)
    dims_MB = (NP_star..., length(molsys.properties.bond_types))

    weight_bond_hat = zeros(Float64, dims_MB)
    weight_bond = zeros(Float64, NP_star)  # Workspace for each bond

    σ = molsys.properties.monomers[:diameters]

    for (u, (i, j)) in enumerate(molsys.properties.bond_types)
        σ_ij = (σ[i] + σ[j]) / 2.0  # Bond center separation
        delta_func_fft!(weight_bond, bin_width, σ_ij)

        for K in CartesianIndices(NP_star)
            weight_bond_hat[K, u] = weight_bond[K] / (4π * σ_ij^2)
        end
    end

    return weight_bond_hat
end

function process_propagator(molsys, geometry::CoordSystem)
    error("No process_propagator method defined for geometry of type $(typeof(geometry))")
end
