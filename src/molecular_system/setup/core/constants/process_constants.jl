"""
    process_constants(dataset::Dict) -> ConstantsStruct

Extract physical constants from the input `dataset`. If a reduced unit length `l_sc` is specified
under `reduced_units`, it is converted to meters. Defaults to 0.5 nm if not provided.

Returns a `ConstantsStruct` containing ε₀, e₀, k_B, N_A, and l_sc.
"""

function process_constants(dataset :: Dict)
    l_sc = haskey(dataset, "reduced_units") && haskey(dataset["reduced_units"], "l_sc") ?
           dataset["reduced_units"]["l_sc"] * 1e-9 :
           0.5e-9

    return ConstantsStruct(l_sc = l_sc)
end

function process_constants(dataset)
    throw(ArgumentError("process_constants expects a Dict, but got $(typeof(dataset))"))
end
