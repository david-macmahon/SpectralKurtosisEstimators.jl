module SpectralKurtosisEstimators

using SpecialFunctions
using Distributions
using ApproxFun
using Memoize
using Roots

import Statistics: mean, std, var, quantile
import StatsBase: skewness, kurtosis
import Distributions: cdf, pdf

export SKEstimator, skhat, s1s2, s1s2!, pearson_criterion, pearson_distribution
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

# Size of the s1/s2/sk output arrays for reductions of A along dims over
# blocks of M*N samples; also errors if there are too few samples
function _s1s2_outsize(A::AbstractArray, M, N, dims)
    T = length(axes(A, dims)) ÷ (M*N)
    T > 0 || error("$(size(A, dims)) is too few samples for M=$(M) and N=$(N)")
    (size(A)[1:dims-1]..., T, size(A)[dims+1:end]...), T
end

"""
    s1s2(A, M, N; dims=ndims(A))
    s1s2(A, ske::SKEstimator; dims=ndims(A))
    s1s2(A; dims=ndims(A))

Compute the summed power `s1` and the sum of squared accumulated power `s2`
of the power data in `A` along dimension `dims` for `M` outer sum addends and
`N` inner sum addends.  These are the `s1` and `s2` values consumed by
`skhat`, i.e. `s1 = Σᵢ (Σₖ Pᵢₖ)` and `s2 = Σᵢ (Σₖ Pᵢₖ)²` where `k` indexes
the `N` inner addends and `i` the `M` outer addends of each block.

The computations are fused reductions (via `sum!` and `Base.mapreducedim!`)
that never materialize a full-size intermediate array.

Returns named tuple `(; s1, s2)`, where `s1` and `s2` have the same dimensions
as `A` except that dimension `dims` will be `size(A, dims) ÷ (M*N)`.  If
`M*N` does not divide `size(A, dims)` evenly then some number (less than
`M*N`) of samples from the end of the `dims` dimension of `A` will not be
used (as in `skhat`).  Use [`s1s2!`](@ref) to write into preallocated buffers.

If `M` and `N` are omitted, the entire `dims` dimension of `A` is collapsed
into a single accumulation per slice: `M = size(A, dims)` and `N = 1`, so
`s1` and `s2` have a singleton in the `dims` position (`vec`/`dropdims` them
for a true vector).  For blocked estimates, pass `M` explicitly.
"""
function s1s2(A::AbstractArray, M, N; dims::Integer=ndims(A))
    outsize, _ = _s1s2_outsize(A, M, N, dims)
    s1 = similar(A, eltype(A), outsize)
    s2 = similar(A, eltype(A), outsize)
    s1s2!(s1, s2, A, M, N; dims)
end

function s1s2(A::AbstractArray, ske::SKEstimator; dims::Integer=ndims(A))
    s1s2(A, ske.M, ske.N; dims)
end

function s1s2(A::AbstractArray; dims::Integer=ndims(A))
    s1s2(A, size(A, dims), 1; dims)
end

"""
    s1s2!(s1, s2, A, M, N; dims=ndims(A))
    s1s2!(s1, s2, A, ske::SKEstimator; dims=ndims(A))
    s1s2!(s1, s2, A; dims=ndims(A))

Like [`s1s2`](@ref), but writes the results into the preallocated arrays `s1`
and `s2`, overwriting their contents, so that buffers can be reused across
calls for allocation-free processing.  Combined with broadcasted
`skhat.(s1, s2, ske)` into a preallocated `sk` array, this allows fully
allocation-free computation of spectral kurtosis estimates.

`s1` and `s2` must have the same dimensions as `A` except that dimension
`dims` will be `size(A, dims) ÷ (M*N)`; a `DimensionMismatch` is thrown
otherwise.  Returns `(; s1, s2)`.

If `M` and `N` are omitted, the entire `dims` dimension of `A` is collapsed
into a single accumulation per slice (`M = size(A, dims)` and `N = 1`), and
`s1` and `s2` must have a singleton in the `dims` position.
"""
function s1s2!(s1, s2, A::AbstractArray, M, N; dims::Integer=ndims(A))
    M = Int(M)
    N = Int(N)
    M >= 1 || throw(ArgumentError("M ($M) must be >= 1"))
    N >= 1 || throw(ArgumentError("N ($N) must be >= 1"))
    outsize, T = _s1s2_outsize(A, M, N, dims)
    size(s1) == outsize ||
        throw(DimensionMismatch("size(s1) = $(size(s1)) must be $outsize"))
    size(s2) == outsize ||
        throw(DimensionMismatch("size(s2) = $(size(s2)) must be $outsize"))

    axesA = axes(A)
    Ausable = if M*N*T == size(A, dims)
        A
    else
        view(A, axesA[1:dims-1]..., 1:M*N*T, axesA[dims+1:end]...)
    end

    pre = size(A)[1:dims-1]
    post = size(A)[dims+1:end]
    Anmt = reshape(Ausable, pre..., N, M, T, post...)

    # Inner sums over the N addends: for N == 1 the reshaped input already has
    # the M/T dims in place, so no pass over the data is needed.  The result is
    # kept in an (M, T, ...) layout with all singleton dims dropped so that the
    # outer reductions below hit exactly one dimension each, matching the
    # reduction order of sum(...; dims)/sum(abs2, ...; dims) bit for bit
    if N == 1
        s0 = reshape(Anmt, pre..., M, T, post...)
    else
        s0 = similar(Anmt, pre..., M, T, post...)
        sum!(reshape(s0, pre..., 1, M, T, post...), Anmt)
    end

    # s1: sum of inner sums over the M outer addends
    sum!(reshape(s1, pre..., 1, T, post...), s0)
    # s2: sum of squared inner sums over the M outer addends
    # (Base.mapreducedim! is the internal workhorse behind sum(x -> f(x), A; dims)
    # used here so that the result lands in the caller-provided s2 buffer;
    # it accumulates into s2 without initializing it, so zero s2 first)
    fill!(s2, zero(eltype(s2)))
    Base.mapreducedim!(x -> x^2, +, reshape(s2, pre..., 1, T, post...), s0)

    (; s1, s2)
end

function s1s2!(s1, s2, A::AbstractArray, ske::SKEstimator; dims::Integer=ndims(A))
    s1s2!(s1, s2, A, ske.M, ske.N; dims)
end

function s1s2!(s1, s2, A::AbstractArray; dims::Integer=ndims(A))
    s1s2!(s1, s2, A, size(A, dims), 1; dims)
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

This method allocates the output Arrays each call.  For allocation-free
processing, use [`s1s2!`](@ref) to write `s1`/`s2` into preallocated buffers
and broadcast `skhat.(s1, s2, ske)` into a preallocated `sk` array.
"""
function skhat(A::AbstractArray, ske::SKEstimator; dims::Integer=ndims(A))
    (; s1, s2) = s1s2(A, ske; dims)
    sk = skhat(s1, s2, ske)

    (; s1, sk)
end

end # module SpectralKurtosisEstimators
