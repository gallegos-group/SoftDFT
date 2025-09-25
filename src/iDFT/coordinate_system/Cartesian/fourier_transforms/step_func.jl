function step_func_fft!(weight_fft, DR, bval)
    NP = size(weight_fft)
    dims = ndims(weight_fft)
    Rsys = CartesianIndices(NP)

    if dims != length(DR)
        error("The number of dimensions of weight_fft must match the length of DR")
    end

    @. weight_fft = 0.0
    weight_fft[first(Rsys)] = (4.0 / 3.0) * pi * bval^3

    for K in Iterators.drop(Rsys, 1)
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
        weight_fft[K] = 4.0 * pi * (sin(bval * dist) - bval * dist * cos(bval * dist)) / dist^3
    end

    return weight_fft
end
