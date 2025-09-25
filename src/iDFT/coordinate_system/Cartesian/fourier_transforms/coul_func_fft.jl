function coul_func_fft!(weight_fft, DR, BJ)
    NP = size(weight_fft)
    dims = ndims(weight_fft)
    Rsys = CartesianIndices(NP)

    @. weight_fft = 0.0  # Zero out all entries including the origin

    for K in Iterators.drop(Rsys, 1)  # Skip the zero-frequency point
        idx = Tuple(K)
        temp = 0.0

        for v = 1:dims
            jv = idx[v] - 1
            N = NP[v]

            if jv <= div(N, 2)
                kv = 2.0 * pi * jv / (N * DR[v])
            else
                kv = -2.0 * pi * (N - jv) / (N * DR[v])
            end

            temp += (1.0 - cos(kv * DR[v])) / DR[v]^2
        end

        weight_fft[K] = 2.0 * pi * BJ / temp
    end

    return weight_fft
end

