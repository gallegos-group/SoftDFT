using YAML
using DelimitedFiles
using Printf

include("src/iDFT/iDFT.jl")

# Helper: consistent folder names like dist_1.60, dist_1.62, ...
_format_dir(d; digits=2) = "pH_" * @sprintf("%0.*f", digits, d)

function adjust_pH(input::String; startZ=3.0, stopZ=10.0, stepZ=0.05, out_root::String="pH_sweep")
    data = YAML.load_file(input)
    return adjust_pH(data; startZ=startZ, stopZ=stopZ, stepZ=stepZ, out_root=out_root)
end
# assumes: _format_dir, read_brush_height, read_avg_state_fraction already defined

function adjust_pH(data; startZ=3.0, stopZ=10.0, stepZ=0.05, out_root::String="pH_sweep")
    results = NamedTuple[]  # (step, pH, brush_height, alpha)

    root_abs = abspath(out_root)  # <-- absolute root
    mkpath(root_abs)
    original_pwd = pwd()

    old_dft_system = nothing
    dft_system     = nothing

    try
        for (i, pH) in enumerate(startZ:stepZ:stopZ)
            pH_dir = joinpath(out_root, _format_dir(pH; digits=2))
            mkpath(pH_dir)

            cd(pH_dir) do
                # update pSpecies for "Y"
                idxY = findfirst(sp -> get(sp, "sequence", nothing) == "Y", data["species"])
                if idxY !== nothing
                    data["species"][idxY]["pSpecies"] = pH
                end

                # build systems
                bulk_system = iLST(data)
                geometry    = geometry_system(data)
                dft_system  = setup_iDFT(data, bulk_system, geometry)

                if i == 1
                    eval_inhomogeneous(dft_system)
                else
                    process_guess(dft_system)

                    NP       = dft_system.geometry.NP
                    mirrored = dft_system.geometry.features[:mirrored]
                    NP_star  = compute_half_domain(NP, mirrored)
                    Rsys     = CartesianIndices(NP_star)

                    rho_segments_K = dft_system.fields.rho_K.segments
                    rho_bonds_K    = dft_system.fields.rho_K.bonds

                    keep_segments  = old_dft_system.fields.rho_K.segments
                    keep_bonds     = old_dft_system.fields.rho_K.bonds

                    for K in Rsys
                        tK = Tuple(K)
                        if tK[end] < round(Int, NP_star[end]/2) + 2
                            K1 = CartesianIndex(tK[1:end-1]..., NP_star[end] - tK[end] + 1)

                            for u in eachindex(keep_segments)
                                @. rho_segments_K[u][K,   :, :] = keep_segments[u][K,   :, :]
                                @. rho_segments_K[u][K1,  :, :] = keep_segments[u][K,   :, :]
                            end

                            if size(rho_bonds_K, 3) > 0
                                rho_bonds_K[K,  :, :] = keep_bonds[K,  :, :]
                                rho_bonds_K[K1, :, :] = keep_bonds[K,  :, :]
                            end
                        end
                    end

                    @. dft_system.fields.excess.PsiC = old_dft_system.fields.excess.PsiC
                    for rho_segments_u in rho_segments_K
                        apply_mirroring!(rho_segments_u, geometry)
                    end
                    apply_mirroring!(rho_bonds_K, geometry)

                    eval_inhomogeneous(dft_system, 1)
                end
                
                # --- parse outputs ---
                brush_by_species = read_brush_height("brush_height.txt")
                brush_height_val = isempty(brush_by_species) ? NaN : only(values(brush_by_species))

                α = read_avg_state_fraction("total_state_fractions_AAAAAAAAAAAAAAAAAAAAAAAAA.txt"; state_index=2)
                alpha_val = (α === nothing) ? NaN : α

                r = (step=i, pH=pH, brush_height=brush_height_val, alpha=alpha_val)
                push!(results, r)

                _append_summary_line(root_abs, r)  # <-- pass absolute root

                # per-folder info (unchanged)
                mode = (i == 1) ? "w" : "a"
                open("run_info.csv", mode) do io
                    if i == 1
                        println(io, "step,pH,brush_height,alpha")
                    end
                    println(io, "$i,$pH,$brush_height_val,$alpha_val")
                end

                println("Step $i | pH = $pH | α = $alpha_val")

                if alpha_val ≥ 0.99 && sign(stepZ) > 0.0
                    println("Stopping early: α reached $(alpha_val) at step $i (pH = $pH)")
                    return results
                elseif alpha_val ≤ 0.01 && sign(stepZ) < 0.0
                    println("Stopping early: α reached $(alpha_val) at step $i (pH = $pH)")
                    return results
                end
            end

            old_dft_system = deepcopy(dft_system)
        end
    finally
        cd(original_pwd)
    end

    return results
end


# function adjust_pH(data; startZ=3.0, stopZ=10.0, stepZ=0.05, out_root::String="pH_sweep")
#     results = NamedTuple[]  # (step, pH, brush_height, alpha)

#     mkpath(out_root)
#     original_pwd = pwd()

#     old_dft_system = nothing
#     dft_system     = nothing

#     try
#         for (i, pH) in enumerate(startZ:stepZ:stopZ)
#             pH_dir = joinpath(out_root, _format_dir(pH; digits=2))
#             mkpath(pH_dir)

#             cd(pH_dir) do
#                 # update pSpecies for "Y"
#                 idxY = findfirst(sp -> get(sp, "sequence", nothing) == "Y", data["species"])
#                 if idxY !== nothing
#                     data["species"][idxY]["pSpecies"] = pH
#                 end

#                 # build systems
#                 bulk_system = iLST(data)
#                 geometry    = geometry_system(data)
#                 dft_system  = setup_iDFT(data, bulk_system, geometry)

#                 if i == 1
#                     eval_inhomogeneous(dft_system)
#                 else
#                     process_guess(dft_system)

#                     NP       = dft_system.geometry.NP
#                     mirrored = dft_system.geometry.features[:mirrored]
#                     NP_star  = compute_half_domain(NP, mirrored)
#                     Rsys     = CartesianIndices(NP_star)

#                     rho_segments_K = dft_system.fields.rho_K.segments
#                     rho_bonds_K    = dft_system.fields.rho_K.bonds

#                     keep_segments  = old_dft_system.fields.rho_K.segments
#                     keep_bonds     = old_dft_system.fields.rho_K.bonds

#                     for K in Rsys
#                         tK = Tuple(K)
#                         if tK[end] < round(Int, NP_star[end]/2) + 2
#                             K1 = CartesianIndex(tK[1:end-1]..., NP_star[end] - tK[end] + 1)

#                             for u in eachindex(keep_segments)
#                                 @. rho_segments_K[u][K,   :, :] = keep_segments[u][K,   :, :]
#                                 @. rho_segments_K[u][K1,  :, :] = keep_segments[u][K,   :, :]
#                             end

#                             if size(rho_bonds_K, 3) > 0
#                                 rho_bonds_K[K,  :, :] = keep_bonds[K,  :, :]
#                                 rho_bonds_K[K1, :, :] = keep_bonds[K,  :, :]
#                             end
#                         end
#                     end

#                     @. dft_system.fields.excess.PsiC = old_dft_system.fields.excess.PsiC
#                     for rho_segments_u in rho_segments_K
#                         apply_mirroring!(rho_segments_u, geometry)
#                     end
#                     apply_mirroring!(rho_bonds_K, geometry)

#                     eval_inhomogeneous(dft_system, 1)
#                 end

#                 # --- parse outputs ---
#                 brush_by_species = read_brush_height("brush_height.txt")
#                 brush_height_val = isempty(brush_by_species) ? NaN : only(values(brush_by_species))

#                 α = read_avg_state_fraction("total_state_fractions_AAAAAAAAAAAAAAAAAAAAAAAAA.txt"; state_index=2)
#                 alpha_val = (α === nothing) ? NaN : α

#                 # single, consistent push:
#                 push!(results, (step=i, pH=pH, brush_height=brush_height_val, alpha=alpha_val))

#                 # per-folder summary (header only once)
#                 mode = (i == 1) ? "w" : "a"
#                 open("run_info.csv", mode) do io
#                     if i == 1
#                         println(io, "step,pH,brush_height,alpha")
#                     end
#                     println(io, "$i,$pH,$brush_height_val,$alpha_val")
#                 end

#                 println("Step $i | pH = $pH")
#             end

#             old_dft_system = deepcopy(dft_system)
#         end
#     finally
#         cd(original_pwd)
#     end

#     # consolidated summary
#     summary_csv = "summary_brush_alpha.csv"
#     open(summary_csv, "w") do io
#         println(io, "step,pH,brush_height,alpha")
#         for r in results
#             println(io, "$(r.step),$(r.pH),$(r.brush_height),$(r.alpha)")
#         end
#     end

#     return results
# end



# --- Parsers ---------------------------------------------------------------

"""
    read_brush_height(path="brush_height.txt")

Returns a Dict{String,Float64} mapping species name -> brush height.
If there is exactly one species line, also returns that scalar as `val`
for convenience: `val = only(values(dict))`.
"""
function read_brush_height(path::String="brush_height.txt")
    if !isfile(path)
        return Dict{String,Float64}()  # no fixed species / no output
    end
    vals = Dict{String,Float64}()
    open(path, "r") do io
        for ln in eachline(io)
            s = strip(ln)
            isempty(s)      && continue
            startswith(s,"#") && continue
            parts = split(s)
            length(parts) < 2 && continue
            species = parts[1]
            h = try
                parse(Float64, parts[2])
            catch
                continue
            end
            vals[species] = h
        end
    end
    return vals
end

"""
    read_avg_state_fraction(path_glob="total_state_fractions_*.txt"; state_index=2)

Scans total_state_fractions files, finds the line starting with "avg",
and returns the chosen state's average fraction (ionization) per file.
Returns a Dict{String,Float64} mapping species -> α (avg fraction for that state).
`state_index` is 1-based on the `avg` line order.
"""
function read_avg_state_fraction(file::String; state_index::Int=2)
    α = nothing
    open(file, "r") do io
        for ln in eachline(io)
            s = strip(ln)
            startswith(s, "avg ") || continue
            parts = split(s)
            if length(parts) >= 1 + state_index
                val_str = parts[1 + state_index]
                if !isempty(val_str)
                    α = parse(Float64, val_str)
                end
            end
        end
    end
    return α
end

# --- Helper: append one line to the run summary (creates file with header if needed)
# root_abs must be absolute
function _append_summary_line(root_abs::AbstractString, r)
    path = joinpath(root_abs, "summary_brush_alpha.csv")
    mkpath(root_abs)  # ensure root exists

    header_needed = !isfile(path) || filesize(path) == 0
    open(path, header_needed ? "w" : "a") do io
        if header_needed
            println(io, "step,pH,brush_height,alpha")
        end
        @printf(io, "%d,%.12g,%.12g,%.12g\n", r.step, r.pH, r.brush_height, r.alpha)
        flush(io)
    end
end

