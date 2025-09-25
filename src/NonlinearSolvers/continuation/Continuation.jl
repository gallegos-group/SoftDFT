# function objective_pac(problem_data::Tuple)

#     base_objective = problem_data[end-1]

#     updated_vars = base_objective(problem_data)

#     continuation_data = problem_data[end]

#     X_k = continuation_data[:X_k]
#     t = continuation_data[:t]
#     ds = continuation_data[:ds]

#     # Flatten updated_vars into a single vector
#     n = length(X_k)
#     X_out = zeros(Float64, n)
#     flatten_solution!(X_out, updated_vars)

#     # Add the PAC constraint residual
#     updated_cont = X_k[end] + (ds - dot(t[1:end-1], X_out[1:end-1] .- X_k[1:end-1])) / t[end]
    
#     return (updated_vars..., updated_cont)
# end



# function  continuation_function(
#                     continuation_type, # :: PAC
#                     solver_type :: AbstractMixing,
#                     base_problem_data :: Tuple,
#                     base_objective :: Function,
#                     base_solution_variables :: Tuple,
#                     continuation_variable,
#                     numerical_details :: Dict = Dict(),
#                     constraints = nothing,
#                     error_file :: String = "error_pac.txt"
#                                 )


#     count_total_entries(var) = 
#     isa(var, Number) ? 1 :
#     isa(var, AbstractArray) && !isempty(var) && isa(first(var), Number) ? length(var) :
#     isa(var, AbstractArray) && isempty(var) ? 0 :
#     isa(var, AbstractVector) ? sum(count_total_entries.(var)) :
#     isa(var, ContinuationVariable) ? 1 :
#     error("Unsupported structure in solution_variables")

#     solution_variables = (base_solution_variables..., continuation_variable)
#     ni = map(count_total_entries, solution_variables)
#     n = sum(ni)

#     # Default constraints are identity (no change)
#     if constraints === nothing
#         constraints = fill(identity, length(solution_variables))
#     end
#     base_constraints = constraints[1:end-1]

#     num_keep = get(numerical_details, "num_keep", 5)
#     ds_init = get(numerical_details, "ds_init", 0.05)
#     ds_min = get(numerical_details, "ds_min", 1e-4)
#     ds_max = get(numerical_details, "ds_max", 1e-1)
#     gamma = get(numerical_details, "gamma", 0.1)
#     n_target = get(numerical_details, "n_target", 20)
#     max_steps = get(numerical_details, "max_steps", 20000)
#     print_frequency = get(numerical_details, "print_frequency", 1)
#     tole = get(numerical_details, "tole", 1e-6)
#     target_continuation = get(numerical_details, "target_continuation", nothing)

#     # === Storage and Initialization ===
#     X_k = zeros(n)
#     X_prev = zeros(n)
#     X_pred = zeros(n)
#     X_hist = zeros(n, num_keep)
#     S_vec = zeros(num_keep)
#     dYdS = zeros(n)
#     dS = ds_init
#     s = 0.0

#     # Solve for initial continuation value
#     starting_point = get_value(solution_variables[end])
#     println("\nInitial solution for $starting_point")
#     solver_function(
#         solver_type, 
#         base_problem_data,
#         base_objective,
#         base_solution_variables,
#         numerical_details,
#         base_constraints
#         )

#     flatten_solution!(X_k, solution_variables)

#     @. X_prev = X_k

#     # Solve for the next continuation value
#     increment!(solution_variables[end], -dS) # small initial step

#     println("\nSecond solution for $(get_value(solution_variables[end]))")
    
#     solver_function(
#         solver_type, 
#         base_problem_data,
#         base_objective,
#         base_solution_variables,
#         numerical_details,
#         base_constraints
#         )

#     flatten_solution!(X_k, solution_variables)

#     # Now update s by the actual distance
#     ds_actual = norm(X_k .- X_prev)
#     s += ds_actual

#     @. X_hist[:, 1] = X_prev
#     @. X_hist[:, 2] = X_k

#     println(X_hist[end,1], "   ", X_hist[end, 2])
#     S_vec[1] = 0.0
#     S_vec[2] = s

#     continuation_data = Dict(:X_k => X_k, :t => dYdS, :ds => dS)
#     problem_data = (base_problem_data..., base_objective, continuation_data)
#     objective_function = objective_pac


#     step = 3
#     # Open output file for writing
#     open("pac_output.txt", "w") do io

#         while pursuing_target(get_value(solution_variables[end]), starting_point, target_continuation;
#                             inclusive=true, check_direction=true) &&
#                             step <= max_steps

#             i1, i2, i3 = mod1(step - 1, num_keep), mod1(step - 2, num_keep), mod1(step - 3, num_keep)

#             if step == 3 || num_keep < 3
#                 @. dYdS = X_hist[:, i1] - X_hist[:, i2]
#             else
#                 ds1 = S_vec[i1] - S_vec[i2]
#                 ds2 = S_vec[i2] - S_vec[i3]
#                 a = ds1 / (ds2 * (ds1 + ds2))
#                 b = -(ds1 + ds2) / (ds1 * ds2)
#                 c = (2.0 * ds1 + ds2) / (ds1 * (ds1 + ds2))
#                 @. dYdS = a * X_hist[:, i3] + b * X_hist[:, i2] + c * X_hist[:, i1]
#             end

#             nrm = norm(dYdS)
#             @. dYdS /= nrm

#             # Save previous converged solution
#             @. X_prev = X_k

#             # Predict next solution
#             @. X_pred = X_k + dS * dYdS
#             apply_flat_solution!(solution_variables, X_pred)

#             println("\n$step solution for $(get_value(solution_variables[end]))")

#             # Calculate alpha (specific property from your solution)
#             alpha = sum(solution_variables[3][1][1, :]) / sum(solution_variables[3][1])
#             height = get_height_brush(problem_data)

#             # Solve the system starting from predicted solution
#             n_iter = solver_function(
#                             solver_type, 
#                             problem_data,
#                             objective_function,
#                             solution_variables,
#                             numerical_details,
#                             constraints
#                         )

#             # Flatten updated solution
#             flatten_solution!(X_k, solution_variables)

#             # Compute true step taken
#             ds_actual = norm(X_k .- X_prev)

#             # Update total arclength based on actual step
#             s += ds_actual

#             dS = (dS + ds_actual)/2.0
#             # === Adaptive step size control based on solver iterations ===
#             # factor = (n_target / n_iter)^gamma
#             # dS_new = dS * factor
#             # dS = clamp(dS_new, ds_min, ds_max)

#             # if ds_actual < 0.5 * dS
#             #     dS = max(ds_min, dS * 0.5)
#             # elseif ds_actual > 1.2 * dS
#             #     dS = min(ds_max, dS * 1.1)
#             # end
#             # dS = clamp(dS, ds_min, ds_max)

#             # Store updated solution in history
#             i = mod1(step, num_keep)
#             @. X_hist[:, i] = X_k
#             S_vec[i] = s

#             # Print to screen
#             if step % print_frequency == 0
#                 println("PAC Step $step | Continuation param = $(get_value(solution_variables[end])) | alpha = $alpha | height = $height")
#             end

#             # === NEW: Write output to file ===
#             continuation_param = get_value(solution_variables[end])
#             write(io, "$(continuation_param) $(alpha) $(height)\n")
#             flush(io)
#             step += 1
#         end

#     end  # closes output file after while-loop finishes
# end



# function call_continuation(system_input, continuation_variable)

#     struct_iLST = system_input.struct_iLST

#     geometry = system_input.geometry

#     numerical_details = system_input.numerical_details

#     fields = system_input.fields

#     struct_functionals = system_input.struct_functionals

#     base_problem_data = (struct_iLST, geometry, fields, struct_functionals, )

#     base_objective = objective_dft

#     base_solution_variables = (fields.rho_K.beads, fields.rho_K.bonds, fields.rho_K.states, fields.PsiC)
    
#     constraints = (x -> max(x, 0), x -> max(x, 0),  x -> max(x, 0), identity, identity)

#     println(continuation_variable)
#     return continuation_function(
#                 1.0, # :: PAC
#                 AndersonMixing(),
#                 base_problem_data,
#                 base_objective,
#                 base_solution_variables,
#                 continuation_variable,
#                 numerical_details,
#                 constraints
#                                     )

# end


# get_pH = () -> result.struct_iLST.struct_molsys.properties.system["pH"]
# set_pH = (val) -> (result.struct_iLST.struct_molsys.properties.system["pH"] = val)

# struct ContinuationVariable
#     getter::Function
#     setter::Function
#     updater::Function  # <-- new: a post-update hook
# end

# function increment!(cv::ContinuationVariable, delta)
#     set_value!(cv, get_value(cv) + delta)
# end

# function get_value(cv::ContinuationVariable)
#     return cv.getter()
# end

# function set_value!(cv::ContinuationVariable, new_value)
#     cv.setter(new_value)
#     cv.updater()  # <-- run post-update recomputation
# end

# function update_pH_dependent_terms1(result)
#     properties_dict = result.struct_iLST.struct_molsys.properties

#     for (i, isrxn) in enumerate(properties_dict.reactions["participating"])
#         if isrxn
#             pK = properties_dict.monomers["pK_value"][i]
#             pH = properties_dict.system["pH"]
#             properties_dict.monomers["mu_protonation"][i] = (pK - pH) * log(10.0)
#         end
#     end
# end

# pH_continuation = ContinuationVariable(get_pH, set_pH, () -> update_pH_dependent_terms1(result))