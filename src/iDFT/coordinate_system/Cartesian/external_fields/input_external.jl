using DelimitedFiles

#   external_field:
    # input_external

function eval_External(molsys, geometry::CartesianCoord, fields, ::Val{:input_external})
    @unpack Ext, trapez = fields.excess
    @unpack NP = geometry
    NP = geometry.NP

    n_species = size(Ext, ndims(Ext))  # Ext assumed to be (NP..., n_species)

    description = geometry.features["external_field"]["input_external"]
    filepath = get(description, "file", "external_potential.csv")

    data = readdlm(filepath, ',', Float64)
    n_dim = length(NP)
    expected_cols = n_dim + n_species

    @assert size(data, 2) == expected_cols "Expected $expected_cols columns: dims + $n_species species."

    for row in eachrow(data)
        dims = Tuple(round.(Int, row[1:n_dim]) .+ 1)
        idx = CartesianIndex(dims)
        for s = 1:n_species
            Ext[idx, s] = row[n_dim + s]
        end
    end
end