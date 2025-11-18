"""
    excess_free_energy(molsys, bulk) -> Float64

Computes the total excess Helmholtz free energy by summing contributions
from all model terms except the ideal gas term.

Returns:
- `fe_ex`: the total excess free energy
"""
function excess_free_energy(molsys::MolecularSystem, bulk::BulkState)
    fe_ex = 0.0
    for model_term in bulk.fe_model
        if !is_ideal(model_term)
            fe_ex += free_energy(molsys, bulk, model_term) :: Float64
        end
    end
    return fe_ex
end

function free_energy(molsys::MolecularSystem, bulk::BulkState, model_term)
    error("$model_term has no defined free_energy(molsys, bulk, model_term)")
end

