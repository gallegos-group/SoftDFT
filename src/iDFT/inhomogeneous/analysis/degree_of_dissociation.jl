function write_total_state_fractions(rho_segments::Vector{<:AbstractArray{Float64}},
                                     molsys::MolecularSystem,
                                     NP::NTuple,
                                     out_prefix::String = "total_state_fractions_")

    @assert length(rho_segments) == length(molsys.configurations)

    nd_spatial   = length(NP)
    configs      = molsys.configurations
    config_names = molsys.properties.species[:species]

    for (i, (rho, config)) in enumerate(zip(rho_segments, configs))
        state_family = config.state_family
        n_segments   = length(state_family)

        # Skip species with no multistate segments
        has_multistate = any(length.(state_family) .> 1)
        if !has_multistate
            continue
        end

        species_name = config_names[i]
        filename     = out_prefix * species_name * ".txt"

        max_states          = maximum(length.(state_family))
        sum_fractions       = zeros(Float64, max_states)
        count_contributions = zeros(Int,     max_states)

        # Prebuild spatial slice tuple (:, :, ...) once
        idxs_spatial = ntuple(_ -> Colon(), nd_spatial)

        open(filename, "w") do file
            # Header
            println(file, "segment ", join(["state$(n)" for n in 1:max_states], " "))

            for seg in 1:n_segments
                n_states = length(state_family[seg])
                if n_states == 1
                    continue  # skip non-dissociating segments
                end

                local_vals = Float64[]
                for local_idx in 1:n_states
                    # Sum density over all spatial bins for this state & segment
                    v = sum(@view rho[idxs_spatial..., local_idx, seg])
                    push!(local_vals, v)
                end

                total  = sum(local_vals)
                normed = total > 0 ? (local_vals ./ total) : zeros(n_states)

                # Update running sums for averaging (state-wise)
                @inbounds for j in 1:n_states
                    sum_fractions[j]       += normed[j]
                    count_contributions[j] += 1
                end

                # Pad to max_states for aligned columns
                padded = vcat(normed, fill("", max_states - n_states))
                println(file, "$seg ", join(padded, " "))
            end

            # Average across dissociating segments (per state index)
            println(file, "\n# Average state fractions across all dissociating segments:")
            avg_fractions = Any[]
            for j in 1:max_states
                if count_contributions[j] > 0
                    push!(avg_fractions, sum_fractions[j] / count_contributions[j])
                else
                    push!(avg_fractions, "")
                end
            end
            println(file, "avg ", join(avg_fractions, " "))
        end
    end

    return nothing
end