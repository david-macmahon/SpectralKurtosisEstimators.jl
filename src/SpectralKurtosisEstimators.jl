module SpectralKurtosisEstimators

using SpecialFunctions, Distributions
using SpecialFunctions: loggammadiv
import Statistics: mean, var
import StatsBase: skewness, kurtosis

export SKEstimator, skhat, pearson_critereon, pearson_type_iii, pearson_type_vi

include("ske.jl")
include("pearson.jl")

"""
    skhat(s1, s2, M, N=1, d=1)
    skhat(s1, s2, ske::SKEstimator)

Compute the spectral kurtosis estimate from:

- `s1`: sum of power
- `s2`: sum of squared power
- `M`, `ske.M`: number of samples summed (e.g. "off-board")
- `N`, `ske.N`: number of samples pre-summed (e.g "on-board")
- `d`, `ske.d`: shape parameter for original voltage data
  - Use `1/2` for real voltages
  - Use `1` for complex voltages

The formulas used here is from equation 8 of:
"Monthly Notices of the Royal Astronomical Society". 406, L60-L64 (2010)
doi:10.1111/j.1745-3933.2010.00882.x
"""
function skhat(s1, s2, M, N=1, d=1)
    @. (M*N*d+1) * (M*s2 / s1^2 - 1) / (M-1)
end

function skhat(s1, s2, ske::SKEstimator)
    skhat(s1, s2, ske.M, ske.N, ske.d)
end

"""
    skhat(A::AbtractArray, ske::SKEstimator; dims=:)

Compute the spectral kurtosis estimate of `A` along `dims`:

- s1: sum of power
- s2: sum of squared power
- M: number of samples summed (e.g. "off-board")
- N: number of samples pre-summed (e.g "on-board")
- d: shape parameter for original voltage data
  - 1/2 for real voltages
  - 1 for complex voltages

The formulas used here is from equation 8 of:
"Monthly Notices of the Royal Astronomical Society". 406, L60-L64 (2010)
doi:10.1111/j.1745-3933.2010.00882.x
"""
function skhat(A::AbstractArray, ske::SKEstimator; dims=:)
    m = prod(size(A)[dims])
    if m != ske.M
        @warn "number of values summed ($m) != SKEstimator.M ($(ske.M))"
    end
    s1 = sum(A; dims)
    s2 = sum(abs2, A; dims)

    skhat(s1, s2, ske)
end

end # module SpectralKurtosisEstimators