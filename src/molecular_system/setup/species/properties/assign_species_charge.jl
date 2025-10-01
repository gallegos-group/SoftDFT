function assign_species_charge!(configurations::Vector{ConfigurationStruct}, properties::PropertiesStruct)
    if !haskey(properties.monomers, :valences)
        # No valence data: skip assignment
        return
    end

    valences = properties.monomers[:valences]  # Vector{Float64}
    n_species = length(configurations)

    # Initialize or overwrite species_charge vector
    species_charge = zeros(Float64, n_species)

    for i in 1:n_species
        sequence = configurations[i].sequence
        for m in sequence
            species_charge[i] += valences[m]
        end
    end

    properties.species[:species_charge] = species_charge
end