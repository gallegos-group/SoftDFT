# === Base case: empty tuple (no functionals) ===
@inline function eval_fe_functional(
    bulk_system::IsingLST,
    geometry::CoordSystem,
    fields::SpatialFields,
    ::Tuple{}
)
    return 0.0
end

# === Recursive case: process the first functional, then recurse ===
@inline function eval_fe_functional(
    bulk_system::IsingLST,
    geometry::CoordSystem,
    fields::SpatialFields,
    functionals::Tuple
)
    return eval_fe_functional(bulk_system, geometry, fields, functionals[1]) +
           eval_fe_functional(bulk_system, geometry, fields, Base.tail(functionals))
end
