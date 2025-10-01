function contact_value_theorem(dft_system :: IsingDFT)
    return contact_value_theorem(dft_system.bulk_system.molsys, dft_system.geometry, dft_system.fields)
end

function contact_value_theorem(molsys, geometry, fields; filename="cvt.txt")

    contr = zeros(Float64, (2,length(geometry.NP)))
    for (name, _) in geometry.features[:external_field]
        contr .+= contact_value_theorem(molsys, geometry, fields, Val{Symbol(name)}())
    end

    open(filename, "w") do io
        println(io, "# Contact Value Theorem results")
        for d in 1:length(geometry.NP)
            println(io, "dim$(d)  $(contr[1,d])  $(contr[2,d])")
        end
    end

    return contr
end