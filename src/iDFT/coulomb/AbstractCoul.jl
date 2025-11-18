abstract type AbstractCoulomb end

include("models/coulomb.jl")
include("neutrality_shift.jl")
include("coulomb_registry.jl")

# Returns true if this functional type is a long-range Coulomb model
is_coulomb_model(::Type) = false

# Override for Coulomb types
is_coulomb_model(::Type{<:AbstractCoulomb}) = true

function construct_coulomb(coulomb :: AbstractCoulomb,
                              molsys::MolecularSystem,
                              bulk::BulkState,
                              geometry::CartesianCoord)
    
    error("No method defined for $(typeof(coulomb)).")
end


function eval_coulomb_energy(bulk_system :: IsingLST, geometry :: CoordSystem, fields :: SpatialFields, coulomb :: AbstractCoulomb)
    error("No method defined for $(typeof(coulomb)).")
end


function solve_coulomb!(bulk_system, geometry, fields, coulomb::AbstractCoulomb)
    error("No method defined for $(typeof(coulomb)).")
end

function solve_poisson_eq(bulk_system, geometry, fields, coulomb :: AbstractCoulomb)
    error("No method defined for $(typeof(coulomb)).")
end