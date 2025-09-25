include("hard_walls.jl")
include("charged_walls.jl")
include("pointcharge_walls.jl")
include("sw_walls.jl")
include("input_external.jl")

function eval_external_field(dft_system)
    return eval_external_field(dft_system.bulk_system.molsys, dft_system.geometry, dft_system.fields)
end

function eval_external_field(molsys, geometry, fields)
    @unpack Ext, trapez = fields.excess
    @unpack surf_hat, plan_forward = fields.fourier

    @. Ext = 0.0
    @. trapez = 1.0
    @. surf_hat = 0.0 + 0im

    for (name, _) in geometry.features[:external_field]
        eval_External(molsys, geometry, fields, Val{Symbol(name)}())
    end

    surf_hat = plan_forward * surf_hat
end


function contact_value_theorem(molsys, geometry, fields, unknown)
    # If unknown
    println("contact_value_theorem not specified for $unknown.")
    return zeros(Float64, size(geometry.NP))
end