"""
    process_fields(rho::RhoStruct, molsys::MolecularSystem, geometry::CartesianCoord) → SpatialFields

Constructs and returns the full set of field data required for inhomogeneous DFT calculations.
This includes:

- Density fields (beads, bonds, simulated)
- Excess fields (external potentials, electrostatics, non-ideality)
- Fourier transform plans and arrays
- Fixed segment structures for anchored or grafted chains

# Arguments
- `rho`: Preprocessed `RhoStruct` with segment and bond metadata
- `molsys`: Full molecular system (used to determine configurations and fixed species)
- `geometry`: The spatial grid and coordinate system (e.g., Cartesian)

# Returns
- `SpatialFields`: A structure bundling all per-grid and per-segment fields
"""
function process_fields(molsys::MolecularSystem, geometry::CoordSystem)
    NP = geometry.NP
    NB = length(molsys.properties.species.monomers)
    MB = length(molsys.properties.bond_types)

    dims_NB = (NP..., NB)
    dims_MB = (NP..., 2, MB)

    # === Density Fields ===
    rho_beads_K = zeros(Float64, dims_NB)
    rho_pairs_K = zeros(Float64, dims_MB)
    rho_sim_K   = zeros(Float64, dims_NB)
    rho_segments_K = [
            zeros(Float64, 
            (NP..., maximum(length, config.state_family), length(config.sequence))) 
            for config in molsys.configurations]

    density_fields = DensityFields(rho_beads_K, rho_pairs_K, rho_sim_K, 
                                    rho_segments_K)

    # === Excess Fields ===
    mu_ex_K = zeros(Float64, dims_NB)
    lng_K   = zeros(Float64, dims_MB)
    Ext     = zeros(Float64, dims_NB)
    trapez  = ones(Float64, dims_NB)
    Psi     = zeros(Float64, NP)
    PsiC    = [0.0]

    excess_fields = ExcessFields(mu_ex_K, lng_K, Ext, trapez, Psi, PsiC)

    # === FFT Cache ===
    fourier = process_fourier(molsys, geometry)

    # === Fixed Species ===
    fixed = [FixedSpecies(length(cfg.sequence), length(NP)) for cfg in molsys.configurations]

    # === Bundle everything into SpatialFields ===
    return SpatialFields(density_fields, excess_fields, fourier, fixed)
end