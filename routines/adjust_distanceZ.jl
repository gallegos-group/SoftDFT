using YAML
using DelimitedFiles
using Printf

include("src/iDFT/iDFT.jl")

# Helper: consistent folder names like dist_1.60, dist_1.62, ...
_format_dist(d; digits=2) = "dist_" * @sprintf("%0.*f", digits, d)

function adjust_distancesZ(input::String; startZ=1.6, stopZ=20.0, stepZ=0.02, errP=1e-6, out_root::String="dist_sweep")
    data = YAML.load_file(input)
    return adjust_distancesZ(data; startZ=startZ, stopZ=stopZ, stepZ=stepZ, errP=errP, out_root=out_root)
end

function adjust_distancesZ(data; startZ=1.6, stopZ=20.0, stepZ=0.02, errP=1e-6, out_root::String="dist_sweep")
    # Storage for overall sweep
    result  = Float64[]
    results = NamedTuple[]  # (step, distance, Pressure, deltaP, energy)

    mkpath(out_root)  # parent sweep folder

    old_dft_system = nothing
    dft_system     = nothing

    err       = errP
    max_count = 20
    distance  = startZ
    count     = 0
    i         = 0

    # Keep a copy to restore PWD even if something throws
    original_pwd = pwd()
    try
        while count < max_count
            distance += stepZ
            i += 1

            # Set up folder for this distance and run inside it
            dist_dir = joinpath(out_root, _format_dist(distance; digits=2))
            mkpath(dist_dir)

            cd(dist_dir) do
                # Update geometry
                data["geometry"]["dimensions"] = [distance]

                # Optionally save effective input used for this run
                # (Uncomment if your YAML lib supports it)
                # open("effective_input.yaml", "w") do io
                #     YAML.write(io, data)
                # end

                # --- Build systems ---
                bulk_system = iLST(data)
                geometry    = geometry_system(data)
                dft_system  = setup_iDFT(data, bulk_system, geometry)

                if i == 1
                    eval_inhomogeneous(dft_system)
                    rho_segments_K = dft_system.fields.rho_K.segments
                    rho_bonds_K    = dft_system.fields.rho_K.bonds
                else
                    process_guess(dft_system)

                    NP       = dft_system.geometry.NP
                    Rsys     = CartesianIndices(NP)

                    rho_segments_K = dft_system.fields.rho_K.segments
                    rho_bonds_K    = dft_system.fields.rho_K.bonds

                    keep_segments  = old_dft_system.fields.rho_K.segments
                    keep_bonds     = old_dft_system.fields.rho_K.bonds

                    for K in Rsys
                        idx = Tuple(K)
                        if idx[end] < round(Int, NP[end]/2) + 2
                            idx1 = CartesianIndex(idx[1:end-1]..., NP[end] - idx[end] + 1)

                            for u in eachindex(keep_segments)
                                @. rho_segments_K[u][K, :, :]   = keep_segments[u][K, :, :]
                                @. rho_segments_K[u][idx1, :, :] = keep_segments[u][K, :, :]
                            end

                            if size(rho_bonds_K, 3) > 0
                                rho_bonds_K[K, :, :]   = keep_bonds[K, :, :]
                                rho_bonds_K[idx1, :, :] = keep_bonds[K, :, :]
                            end
                        end
                    end

                    @. dft_system.fields.excess.PsiC = old_dft_system.fields.excess.PsiC

                    eval_inhomogeneous(dft_system, 1)
                end

                # --- Measurements ---
                dimension = dft_system.geometry.dimensions[1]
                omega     = grand_potential(dft_system)
                energy    = omega + dft_system.bulk_system.bulk.Pressure[1] * dimension

                Pressure = contact_value_theorem(dft_system)[1]

                if length(results) != 0
                    delta_E = energy - results[end].energy
                else
                    delta_E = 10.0
                end

                deltaP = Pressure - dft_system.bulk_system.bulk.Pressure[1]
                push!(result,  deltaP)
                push!(results, (step=i, distance=distance, Pressure=Pressure, deltaP=deltaP, energy=energy))

                # Per-folder tiny log
                open("run_info.csv", "w") do io
                    println(io, "step,distance,Pressure,deltaP,energy")
                    println(io, "$i,$distance,$Pressure,$deltaP,$energy")
                end

                println("Step $i | Distance = $distance | Pressure = $Pressure | ΔP = $deltaP | Energy = $energy")
                println("Currently on count: $count of $max_count")

                if abs(delta_E) < err
                    count += 1
                else
                    count = 0
                end

                if distance > stopZ
                    println("Trouble converging... ΔP = $(abs(deltaP))")
                    count = max_count + 1
                end
            end # cd(dist_dir)

            # Keep the last system for warm-starting next distance
            old_dft_system = deepcopy(dft_system)
        end

        # --- Write overall summary in parent folder ---
        final_energy = results[end].energy
        final_deltaP = results[end].deltaP
        summary_csv  = "adjust_distancesZ_results.csv"
        open(summary_csv, "w") do io
            println(io, "step,distance,Pressure,deltaP,energy,deltaE")
            for r in results
                deltaE = r.energy - final_energy
                deltaP = r.deltaP - final_deltaP
                println(io, "$(r.step),$(r.distance),$(r.Pressure),$deltaP,$(r.energy),$deltaE")
            end
        end

    finally
        cd(original_pwd)
    end

    return result
end
