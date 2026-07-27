module SpectralKurtosisEstimators

using SpecialFunctions
using Distributions
using ApproxFun
using Memoize
using Roots

import Statistics: mean, std, var, quantile
import StatsBase: skewness, kurtosis
import Distributions: cdf, pdf

export SKEstimator, skhat, pearson_criterion, pearson_distribution
export PearsonTypeVI, PearsonTypeIII, PearsonTypeIV

include("pearson_distributions.jl")
include("ske.jl")

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

The formulas used here are from equation 8 of:
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
    skhat(A::AbstractArray, ske::SKEstimator; dims=ndims(A))

Compute generalized spectral kurtosis estimate of `A` along `dims` as specified
by the `M` and `N` fields of `ske`.  `dims` must be an integer and defaults to
the last dimension of `A`.  `ske.d` should be set to half the number of real
values that were pre-summed into each sample of `A`.  Specifically, `ske.N` is
used here to specify additional "inner-sum" integration, so `d` should include a
factor for any inner-sum addends already included in each sample of `A`.

Returns named tuple `(; s1, sk)`, where `s1` is the summed power and `sk` is the
spectral kurtosis.

`s1` and `sk` will have the same dimensions as `A` except that dimension `dims`
will be `size(A, dims) ÷ (M*N)`.  If `M*N` does not divide `size(A, dims)`
evenly then some number (less than `M*N`) of samples from the end of the `dims`
dimension of `A` will not be used.
"""
function skhat(A::AbstractArray, ske::SKEstimator; dims::Integer=ndims(A))
    axesA = axes(A)
    M = ske.M
    N = ske.N
    T = length(axesA[dims]) ÷ (M*N)
    T > 0 || error("$(size(A,dims)) is too few samples for M=$(M) and N=$(N)")

    Ausable = if M*N*T == size(A, dims)
        A
    else
        view(A, axesA[1:dims-1]..., 1:M*N*T, axesA[dims+1:end]...)
    end

    Anmt = reshape(Ausable, axesA[1:dims-1]..., N, M, T, axesA[dims+1:end]...)

    s0 = dropdims(sum(Anmt; dims); dims)
    s1 = dropdims(sum(s0; dims); dims)
    s2 = dropdims(sum(abs2, s0; dims); dims)

    sk = skhat(s1, s2, ske);

    (; s1, sk)
end

end # module SpectralKurtosisEstimators
