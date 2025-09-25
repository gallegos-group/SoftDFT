include("../../species_model/analyticalEval.jl")

"""
    process_evaluation!(dataset::Dict, bulk::BulkState)

Populates `bulk.evaluation` with evaluation strategies for each species.
Valid evaluation types are:
- `"analytical"`: internally computed (e.g. via DFT)
- `"simulation"`: externally provided or imported

If no `"evaluation"` key is given, defaults to `"analytical"`.
"""
function process_evaluation(dataset::Dict)
    species_data = dataset["species"]
    evaluation = Vector{AbstractEvaluation}(undef, length(species_data))

    for (u, species) in enumerate(species_data)
        eval_type = get(species, "evaluation", "analytical")

        evaluator = eval_type == "analytical"  ? AnalyticalEval() :
                    eval_type == "simulation" ? SimulationEval()  :
                    error("Unknown evaluation method '$eval_type'. Must be 'analytical' or 'simulation'.")

        evaluation[u] = evaluator
    end

    return evaluation
end
