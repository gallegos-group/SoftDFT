function write_scalar_value(value::Real, filename::String; label::Union{Nothing, String}=nothing)
    open(filename, "w") do file
        if label === nothing
            println(file, value)
        else
            println(file, label, ": ", value)
        end
    end
end

function write_field_array(arr::AbstractArray{<:Real}, filename::String, NP::NTuple; names::Union{Nothing, Vector{String}}=nothing)
    open(filename, "w") do file
        spatial_dims = NP
        full_dims = size(arr)
        nd_spatial = length(spatial_dims)

        # Check if it's a scalar field (arr has same shape as NP)
        is_scalar = full_dims == spatial_dims

        # Determine how many values per point (e.g., for vector/multi-field case)
        num_values = is_scalar ? 1 : full_dims[nd_spatial + 1]

        # Header
        if names !== nothing
            index_headers = ["Index_$(i)" for i in 1:nd_spatial]
            println(file, join(index_headers, " "), "  ", join(names, " "))
        end

        # Traverse spatial domain
        for idx in CartesianIndices(spatial_dims)
            coord = Tuple(idx)
            values = is_scalar ? [arr[Tuple(idx)...]] : [arr[Tuple(idx)..., j] for j in 1:num_values]
            println(file, join(coord, " "), "  ", join(values, " "))
        end
    end
end
function write_segment_array(rho_segments::Vector{<:AbstractArray{Float64}}, 
                             molsys::MolecularSystem, 
                             NP::NTuple)

    @assert length(rho_segments) == length(molsys.configurations)

    nd_spatial   = length(NP)
    monomer_names = molsys.properties.species.monomers
    configs       = molsys.configurations
    species_names = molsys.properties.species.species

    for (i, (rho, config)) in enumerate(zip(rho_segments, configs))
        @assert size(rho)[1:nd_spatial] == NP
        n_segments = size(rho, nd_spatial + 2)

        # --- Only output if there are >1 segments
        if n_segments <= 1
            continue
        end

        sequence_str  = species_names[i]
        filename      = "rho_segments_" * sequence_str * ".txt"
        state_family  = config.state_family

        open(filename, "w") do file
            # Header
            index_headers  = ["Index_$(d)" for d in 1:nd_spatial]
            column_headers = String[]
            for seg in 1:n_segments
                for (_, global_state) in enumerate(state_family[seg])
                    push!(column_headers, "seg$(seg)_$(monomer_names[global_state])")
                end
            end
            println(file, join(index_headers, " "), "  ", join(column_headers, " "))

            # Data
            idxs_spatial = ntuple(_ -> Colon(), nd_spatial)
            for idx in CartesianIndices(NP)
                coord  = Tuple(idx)
                values = Float64[]
                for seg in 1:n_segments
                    for local_idx in 1:length(state_family[seg])
                        push!(values, @view(rho[idxs_spatial..., local_idx, seg]) |> x->x[idx])
                    end
                end
                println(file, join(coord, " "), "  ", join(values, " "))
            end
        end
    end

    return nothing
end


function write_state_fraction_profiles(rho_segments::Vector{<:AbstractArray{Float64}},
                                       molsys::MolecularSystem,
                                       NP::NTuple,
                                       out_prefix::String = "state_fraction_")
    @assert length(rho_segments) == length(molsys.configurations)

    nd_spatial     = length(NP)
    monomer_names  = molsys.properties.species.monomers
    configs        = molsys.configurations
    config_names   = molsys.properties.species.species

    for (i, (rho, config)) in enumerate(zip(rho_segments, configs))
        state_family = config.state_family
        n_segments   = length(state_family)

        # --- NEW: skip if no segment has multiple states
        has_multistate = any(length.(state_family) .> 1)
        if !has_multistate
            continue
        end

        # Create output filename
        filename = out_prefix * config_names[i] * ".txt"

        open(filename, "w") do file
            # Header
            index_headers  = ["Index_$(d)" for d in 1:nd_spatial]
            column_headers = String[]
            for seg in 1:n_segments
                for state in state_family[seg]
                    push!(column_headers, "seg$(seg)_$(monomer_names[state])")
                end
            end
            println(file, join(index_headers, " "), "  ", join(column_headers, " "))

            # Loop over spatial domain
            for idx in CartesianIndices(NP)
                coord  = Tuple(idx)
                values = Float64[]

                for seg in 1:n_segments
                    denom = 0.0
                    local_vals = Float64[]
                    for (local_idx, _) in enumerate(state_family[seg])
                        val = rho[coord..., local_idx, seg]
                        push!(local_vals, val)
                        denom += val
                    end
                    if denom > 0
                        append!(values, (v / denom for v in local_vals))
                    else
                        append!(values, zeros(length(local_vals)))
                    end
                end

                println(file, join(coord, " "), "  ", join(values, " "))
            end
        end
    end

    return nothing
end