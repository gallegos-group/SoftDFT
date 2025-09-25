function delta_func_fft!(weight_fft, DR, bval)
    NP = size(weight_fft)
    dims = ndims(weight_fft)
    Rsys = CartesianIndices(NP)

    @. weight_fft = 0.0  # Clear array
    weight_fft[first(Rsys)] = 4.0 * pi * bval^2

    for K in Iterators.drop(Rsys, 1)  # Skip zero-frequency point
        idx = Tuple(K)
        dist2 = 0.0

        for v = 1:dims
            jv = idx[v] - 1
            N = NP[v]

            if jv <= div(N, 2)
                kv = 2.0 * pi * jv / (N * DR[v])
            else
                kv = -2.0 * pi * (N - jv) / (N * DR[v])
            end

            dist2 += kv^2
        end

        dist = sqrt(dist2)
        weight_fft[K] = 4.0 * pi * bval * sin(bval * dist) / dist
    end

    return weight_fft
end
