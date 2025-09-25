"""
    apply_filter!(f_hat, geometry; frac_cut=0.6, p=2)

Applies a smooth spectral filter to `f_hat` in Fourier space.

Arguments:
- `f_hat`: FFT array to be filtered (in-place).
- `geometry`: system geometry containing `NP`, `bin_width`, `features`.
- `frac_cut`: fraction of the Nyquist wavenumber to use as cutoff (default 0.6).
- `p`: power in the exponential filter, `exp(-(k/k_cut)^p)` (default 2 for Gaussian).

Returns:
- The filtered `f_hat` array (modified in place).
"""
function apply_filter!(f_hat, geometry; frac_cut=0.5, p=4)
    @unpack NP, bin_width, features = geometry
    @unpack mirrored, offset = features

    NP_star = compute_full_domain(NP, mirrored, offset)

    # Smallest Nyquist frequency across dimensions
    k_nyquist = pi / maximum(bin_width)
    k_cut = frac_cut * k_nyquist

    for K_star in CartesianIndices(f_hat)
        k2 = 0.0
        for d in 1:length(NP)
            n = K_star[d] - 1
            N = NP_star[d]
            Δ = bin_width[d]
            if n > div(N,2)
                n -= N
            end
            kd = 2π * n / (N * Δ)
            k2 += kd * kd
        end
        k = sqrt(k2)
        f_hat[K_star] *= exp(-(k / k_cut)^p)
    end

    return f_hat
end
