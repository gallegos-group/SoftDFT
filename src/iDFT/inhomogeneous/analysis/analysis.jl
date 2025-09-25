include("print_array.jl")
include("grand_potential.jl")
include("degree_of_dissociation.jl")
include("brush_height.jl")
include("output_arrays.jl")
include("contact_value_theorem.jl")


function call_analysis(dft_system :: IsingDFT)

    contact_value_theorem(dft_system)
    write_scalar_value(grand_potential(dft_system), "grand_potential.txt")

    brush_height(dft_system)
    output_arrays(dft_system)
end