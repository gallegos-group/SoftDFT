using Printf

function pac_continuation(input_file :: String)

    dataset = YAML.load_file(input_file)

    return pac_continuation(dataset)
end

function pac_continuation(dataset :: Dict{Any, Any})

    # Solve bulk state (includes mol_sys and bulk_state)
    bulk_system = iLST(dataset)

    # Build geometry
    geom = geometry_system(dataset)

    # Build geometry, fields, functionals
    dft_system = setup_iDFT(dataset, bulk_system, geom)

    eval_external_field(dft_system)

    process_guess(dft_system)

    idxY = findfirst(sp -> get(sp, "sequence", nothing) == "Y", dataset["species"])
    if idxY !== nothing
        cont_var = [dataset["species"][idxY]["pSpecies"]]
    else
        error("Oh no!")
    end

    @unpack bulk_system, geometry, fields, functionals, numerics = dft_system

    problem_data = (bulk_system, geometry, fields, functionals, dataset)

    solution_variables = (fields.rho_K.segments, fields.rho_K.bonds, fields.excess.PsiC)

    continuation_function(problem_data, solution_variables, numerics, cont_var)
end

function continuation_function(problem_data :: Tuple, solution_variables :: Tuple, numerics :: Dict, cont_var :: Vector{Float64})

    ni = map(count_total_entries, solution_variables)
    n = sum(ni) + 1 # for continuation variable

    num_keep = get(numerics, "num_keep", 5)
    ds_init = get(numerics, "ds_init", 0.05)
    ds_min = get(numerics, "ds_min", 1e-4)
    ds_max = get(numerics, "ds_max", 1e-1)
    gamma = get(numerics, "gamma", 0.1)
    n_target = get(numerics, "n_target", 20)
    max_steps = get(numerics, "max_steps", 20000)
    print_frequency = get(numerics, "print_frequency", 1)
    tole = get(numerics, "tole", 1e-6)

    # === Storage and Initialization ===
    X_k = zeros(n)
    X_prev = zeros(n)
    X_pred = zeros(n)
    X_hist = zeros(n, num_keep)
    S_vec = zeros(num_keep)
    dYdS = zeros(n)
    dS = ds_init
    s = 0.0

    solver_function(
                    AndersonMixing(), 
                    problem_data,
                    objective_dft,
                    solution_variables,
                    numerics
                    )

    cont_solution = (solution_variables..., cont_var)
    flatten_solution!(X_k, cont_solution)

    @. X_prev = X_k

    # Solve for next continuation value
    @. cont_var += dS
    dataset = problem_data[end]
    update_cont_var(dataset, cont_var)
    (bulk_system, geometry, fields, functionals, ) = establish_dft_system(dataset)
    problem_data = (bulk_system, geometry, fields, functionals, dataset)
    solution_variables = (fields.rho_K.segments, fields.rho_K.bonds, fields.excess.PsiC)
    
    solver_function(
                    AndersonMixing(), 
                    problem_data,
                    objective_dft,
                    solution_variables,
                    numerics
                    )

    cont_solution = (solution_variables..., cont_var)
    
    flatten_solution!(X_k, cont_solution)

    # Now update s by the actual distance
    ds_actual = norm(X_k .- X_prev)
    s += ds_actual

    @. X_hist[:, 1] = X_prev
    @. X_hist[:, 2] = X_k

    S_vec[1] = 0.0
    S_vec[2] = s

    continuation_data = Dict(:X_k => X_k,
                             :t => dYdS,
                             :ds => dS)

    problem_data = (problem_data..., cont_var, continuation_data)

    step = 3

    results = NamedTuple[]  # (step, pH, brush_height, alpha)

    open("pac_output.txt", "w") do io
        out_root = "pH_sweep"
        root_abs = abspath(out_root)  # <-- absolute root
        mkpath(root_abs)
        original_pwd = pwd()

        while cont_var[1] > 1.99 && cont_var[1] < 12.01

            i1 = mod1(step - 1, num_keep)
            i2 = mod1(step - 2, num_keep)
            i3 = mod1(step - 3, num_keep)

            if step == 3 || num_keep < 3
                @. dYdS = X_hist[:, i1] - X_hist[:, i2]
            else
                ds1 = S_vec[i1] - S_vec[i2]
                ds2 = S_vec[i2] - S_vec[i3]
                a = ds1 / (ds2 * (ds1 + ds2))
                b = -(ds1 + ds2) / (ds1 * ds2)
                c = (2.0 * ds1 + ds2) / (ds1 * (ds1 + ds2))
                @. dYdS = a * X_hist[:, i3] + 
                            b * X_hist[:, i2] + 
                            c * X_hist[:, i1]
            end

            nrm = norm(dYdS)
            @. dYdS /= nrm

            # Save previous converged solution
            @. X_prev = X_k

            # Predict next solution
            @. X_pred = X_k + dS * dYdS
            apply_flat_solution!(cont_solution, X_pred)

            println("\n$step solution for $(cont_var[1])")

            # Solve the system starting from predicted solution
            solver_function(
                            AndersonMixing(), 
                            problem_data,
                            objective_pac,
                            cont_solution,
                            numerics
                        )

            # Flatten updated solution
            flatten_solution!(X_k, cont_solution)

            # Compute true step taken
            ds_actual = norm(X_k .- X_prev)

            # Update total arclength based on actual step
            s += ds_actual

            dS = (dS + ds_actual)/2.0
            # === Adaptive step size control based on solver iterations ===
            # factor = (n_target / n_iter)^gamma
            # dS_new = dS * factor
            # dS = clamp(dS_new, ds_min, ds_max)

            # if ds_actual < 0.5 * dS
            #     dS = max(ds_min, dS * 0.5)
            # elseif ds_actual > 1.2 * dS
            #     dS = min(ds_max, dS * 1.1)
            # end
            # dS = clamp(dS, ds_min, ds_max)

            # Store updated solution in history
            i = mod1(step, num_keep)
            @. X_hist[:, i] = X_k
            S_vec[i] = s

            # --- NEW: Save results for this step ---
            pH = cont_var[1]
            pH_dir = joinpath(out_root, _format_dir(pH; digits=3))
            mkpath(pH_dir)

            objective_dft(problem_data)
            dft_system = IsingDFT(problem_data[1:4]..., numerics)
            cd(pH_dir) do
                call_analysis(dft_system)
                # open("results.txt", "w") do f
                #     println(f, r)
                # end

                # --- parse outputs ---
                brush_by_species = read_brush_height("brush_height.txt")
                brush_height_val = isempty(brush_by_species) ? NaN : only(values(brush_by_species))

                filename = "total_state_fractions_" * dataset["species"][1]["sequence"] * ".txt"
                α = read_avg_state_fraction(filename; state_index=2)
                alpha_val = (α === nothing) ? NaN : α

                r = (step=step, pH=pH, alpha=alpha_val, brush_height=brush_height_val)
                push!(results, r)

                _append_summary_line(root_abs, r)  # <-- pass absolute root

                # per-folder info (unchanged)
                mode = (step == 3) ? "w" : "a"
                open("run_info.csv", mode) do io
                    if step == 3
                        println(io, "step,pH,alpha,brush_height")
                    end
                    println(io, "$i,$pH,$alpha_val,$brush_height_val")
                end
                
                if alpha_val ≥ 0.99 && sign(dS) > 0.0
                    println("Stopping early: α reached $(alpha_val) at step $i (pH = $pH)")
                    return results
                elseif alpha_val ≤ 0.01 && sign(dS) < 0.0
                    println("Stopping early: α reached $(alpha_val) at step $i (pH = $pH)")
                    return results
                end
            end
            
            # Print to screen
            if step % print_frequency == 0
                println("PAC Step $step | Continuation param = $(cont_var[1])")
            end


            # # === NEW: Write output to file ===
            # continuation_param = get_value(solution_variables[end])
            # write(io, "$(continuation_param) $(alpha) $(height)\n")
            # flush(io)
            step += 1
        end

    end  # closes output file after while-loop finishes

    return results
end

function objective_pac(problem_data :: Tuple)

    old_fields = problem_data[3]

    dataset = problem_data[end-2]

    cont_var = problem_data[end-1]

    continuation_data = problem_data[end]

    update_cont_var(dataset, cont_var)

    (bulk_system, geometry, fields, functionals, ) = establish_dft_system(dataset)
    
    for (u, rho_segments_u) in enumerate(fields.rho_K.segments)
        @. rho_segments_u = old_fields.rho_K.segments[u]
    end
    @. fields.rho_K.bonds = old_fields.rho_K.bonds
    @. fields.excess.PsiC = old_fields.excess.PsiC

    updated_vars = objective_dft((bulk_system, geometry, fields, functionals, ))

    X_k = continuation_data[:X_k]
    t = continuation_data[:t]
    ds = continuation_data[:ds]

    # Flatten updated_vars into a single vector
    n = length(X_k)
    X_out = zeros(Float64, n-1)
    flatten_solution!(X_out, updated_vars)

    # Add the PAC constraint residual
    updated_cont = X_k[end] + (ds - dot(t[1:end-1], X_out .- X_k[1:end-1])) / t[end]

    return (updated_vars..., [updated_cont])
end

function update_cont_var(dataset, cont_var:: Vector{Float64})
    idxY = findfirst(sp -> get(sp, "sequence", nothing) == "Y", dataset["species"])
    if idxY !== nothing
        dataset["species"][idxY]["pSpecies"] = cont_var[1]
    else
        error("Oh no!")
    end
end

function establish_dft_system(dataset)

    # Solve bulk state (includes mol_sys and bulk_state)
    bulk_system = iLST(dataset)

    # Build geometry
    geometry = geometry_system(dataset)

    # Build geometry, fields, functionals
    new_dft_system = setup_iDFT(dataset, bulk_system, geometry)

    eval_external_field(new_dft_system)

    return bulk_system, geometry, 
            new_dft_system.fields, 
            new_dft_system.functionals
end




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

# Helper: consistent folder names like dist_1.60, dist_1.62, ...
_format_dir(d; digits=2) = "pH_" * @sprintf("%0.*f", digits, d)