abstract type AbstractEvaluation end

struct AnalyticalEval <: AbstractEvaluation end
struct SimulationEval <: AbstractEvaluation end