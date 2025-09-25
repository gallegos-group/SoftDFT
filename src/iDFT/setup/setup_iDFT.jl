
include("setup_steps/process_external.jl")
include("setup_steps/process_numerical.jl")
include("setup_steps/process_functionals.jl")
include("setup_steps/process_fields.jl")
include("setup_steps/process_fixed.jl")
include("setup_steps/process_guess.jl")

function setup_iDFT(dataset, bulk_system, geometry)

    process_external(bulk_system.molsys, geometry)

    # numerical information
    numerics = process_numerical(dataset)

    # functional information
    fe_functionals = process_functionals(bulk_system, geometry)

    # fields information
    fields = process_fields(bulk_system.molsys, geometry)

    # Process fixed species
    process_fixed(dataset, bulk_system.molsys, fields, geometry)

    eval_external_field(bulk_system.molsys, geometry, fields)

    process_guess(bulk_system, geometry, fields)

    return IsingDFT(bulk_system, geometry, fields, fe_functionals, numerics)
end