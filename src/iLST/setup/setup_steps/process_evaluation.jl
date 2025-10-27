include("../../species_model/AbstractEvaluation.jl")

"""
    process_evaluation!(dataset::Dict, bulk::BulkState)

Populates `bulk.evaluation` with evaluation strategies for each species.
Valid evaluation types are:
- `"analytical"`: internally computed (e.g. via DFT)
- `"simulation"`: externally provided or imported

If no `"evaluation"` key is given, defaults to `"analytical"`.
"""
# function process_evaluation(dataset::Dict)
#     species_data = dataset["species"]
#     evaluation = Vector{AbstractEvaluation}()

#     for (u, species) in enumerate(species_data)
#         eval_type = get(species, "evaluation", "analytical")

#         push!(evaluation, eval_type == "analytical"  ? AnalyticalEval() :
#                     eval_type == "simulation" ? SimulationEval()  :
#                     error("Unknown evaluation method '$eval_type'. Must be 'analytical' or 'simulation'."))
#     end

#     return evaluation
# end
function process_evaluation(dataset::Dict)
    species_data = dataset["species"]
    nt = length(species_data)

    return ntuple(i -> begin
        eval_type = get(species_data[i], "evaluation", "analytical")
        if eval_type == "analytical"
            AnalyticalEval()
        elseif eval_type == "simulation"
            SimulationEval()
        else
            error("Unknown evaluation method '$eval_type'. Must be 'analytical' or 'simulation'.")
        end
    end, nt)
end
