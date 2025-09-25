
include("BulkDensities.jl")
include("../free_energy/Abstract_Types/AbstractFreeEnergy.jl")
include("../species_model/AbstractEvaluation.jl")

"""
    BulkState(rho, struct_model; Psi = 0.0)

Container for the bulk thermodynamic problem and its solution state, including
densities, free energy model components, and thermodynamic output fields.

# Arguments
- `rho`: A `BulkDensities` containing species, segment, bead, and bond densities.
- `struct_model`: A vector of free energy model objects (subtypes of `AbstractFreeEnergy`).
- `Psi`: Optional initial value for the electrostatic potential (either a scalar or a length-1 vector, default = 0.0).

# Fields
- `rho`: The input `BulkDensities` object defining the density fields.
- `mu_species`: Vector of chemical potentials for each species.
- `mu_segments`: Vector of matrices; each matrix corresponds to a species and stores chemical potentials 
  per segment and per internal state. Matches the shape of `rho.segments`.
- `mu_ex`: Vector of excess chemical potentials for each bead (monomer type).
- `lng`: Vector of logarithmic contact values of the radial distribution function, one per bond type.
- `Pressure`: Singleton vector storing the scalar pressure of the system.
- `Psi`: Singleton vector storing the electrostatic potential.
- `struct_model`: Vector of free energy model components used to compute thermodynamic quantities.
"""

struct BulkState
    rho          :: BulkDensities
    mu_species   :: Vector{Float64}
    mu_ex        :: Vector{Float64}
    lng          :: Vector{Float64}
    Xi           :: Vector{Float64}
    Pressure     :: Vector{Float64}
    Psi          :: Vector{Float64}
    fe_model     :: Vector{AbstractFreeEnergy}
    evaluation   :: Vector{AbstractEvaluation}

    function BulkState(rho::BulkDensities, fe_model::Vector{AbstractFreeEnergy}, evaluation :: Vector{AbstractEvaluation}; Psi = 0.0)
        mu_species = zeros(length(rho.species))
        mu_ex      = zeros(length(rho.beads))
        lng        = zeros(size(rho.bonds, 2))
        Xi         = zeros(length(rho.species))
        Pressure   = [0.0]
    
        # Handle Psi input (scalar or singleton vector)
        Psi_vec = isa(Psi, Number) ? [Psi] :
                  (length(Psi) == 1 ? Psi : error("Psi must be a scalar or a length-1 vector."))
    
        return new(rho, mu_species, mu_ex, lng, Xi, Pressure, Psi_vec, fe_model, evaluation)
    end    
end

function Base.show(io::IO, state::BulkState)
    println(io, "BulkState:")
    println(io, "  Species:  ", length(state.rho.species))
    println(io, "  Beads:    ", length(state.rho.beads))
    println(io, "  Pressure: ", state.Pressure[1])
    println(io, "  Psi:      ", state.Psi[1])
    println(io, "  FE Models:   ", join([typeof(m).name.wrapper for m in state.fe_model], ", "))
    println(io, "  Evaluation:   ", join([typeof(m).name.wrapper for m in state.evaluation], ", "))
end