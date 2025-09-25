
model_dir = joinpath(@__DIR__, "species_model")
if !isdir(model_dir)
    error("Required directory 'species_model' not found at $model_dir. Please ensure the models directory exists and contains species_model files.")
end

for file in sort(readdir(model_dir; join=true))
    if endswith(file, ".jl")
        include(file)
    end
end



"""
    determine_species_model(species::Dict, numbeads::Int) -> AbstractSpecies

Given a species dictionary and number of beads in the chain, determines which 
species model to use. If the `species_model` key is specified in the species dictionary, 
that species_model is instantiated with species_model(numbeads). If unspecified, defaults 
to the `freely_jointed` model.

Throws an error if the model name is unrecognized.
"""

function determine_species_model(species, numbeads)
    if haskey(species, "model")
        model_name = species["model"]
        try
            species_model = eval(Symbol(model_name))(numbeads)
        catch
            error("Unknown chain model $(model_name).")
        end
    else
        species_model = freely_jointed(numbeads)
    end

    return species_model
end
