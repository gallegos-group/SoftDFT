function vectored_delta_func_fft!(weight_fft, index::Integer, DR, bval)
    NP = size(weight_fft)
    dims = ndims(weight_fft)
    Rsys = CartesianIndices(NP)

    if dims != length(DR)
        error("The number of dimensions of weight_fft must match the length of DR")
    elseif dims < index
        @. weight_fft = 0.0
        return weight_fft
    else
        weight_fft[Rsys[1]] = 0.0

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

            jv = idx[index] - 1
            N_idx = NP[index]
            if jv <= div(N_idx, 2)
                kv = 2.0 * pi * jv / (N_idx * DR[index])
            else
                kv = -2.0 * pi * (N_idx - jv) / (N_idx * DR[index])
            end

            weight_fft[K] = -4.0 * pi * (sin(bval * dist) - bval * dist * cos(bval * dist)) / dist^3 * kv
        end

        return weight_fft
    end
end
