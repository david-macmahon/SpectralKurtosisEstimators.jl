"""
Base type for all PearsonDistributions
"""
abstract type PearsonDistribution end

# Allow broadcasting over PearsonDistribution objects
Base.broadcastable(d::PearsonDistribution) = Ref(d)

"""
Supertype for types that use approximated functions to compute `pdf`, `cdf`, and
`quantile`.  `PearsonTypeIV` is an (the only) example of this.
"""
abstract type PearsonApproximatedDistribution <: PearsonDistribution end

"""
Supertype for types that can be represented by a known `Distribution`.  For
example, the `PearsonTypeIII` distribution is a location shifted `Gamma`
distribution and the `PearsonTypeVI` distribution is a location shifted
`BetaPrime` distribution.
"""
abstract type PearsonAnalyticDistribution <: PearsonDistribution end

# Concrete types

include("pearson_type_iii.jl")
include("pearson_type_iv.jl")
include("pearson_type_vi.jl")

"""
Compute the Pearson criterion from the second, third, and fourth central
moments `u2`, `u3`, `u4`, resp.
"""
function pearson_criterion(u2, u3, u4)
    B1 = u3^2 / u2^3
    B2 = u4 / u2^2

    (
        (B1 * (B2+3)^2)
        /
        (4(4B2-3B1) * (2B2-3B1-6))
    )
end

# Add statistical/distribution methods

mean(d::PearsonAnalyticDistribution) = mean(distribution(d))
var(d::PearsonAnalyticDistribution) = var(distribution(d))
skewness(d::PearsonAnalyticDistribution) = skewness(distribution(d))
kurtosis(d::PearsonAnalyticDistribution) = kurtosis(distribution(d))

pdf(d::PearsonAnalyticDistribution, x) = pdf(distribution(d), x)
cdf(d::PearsonAnalyticDistribution, x) = cdf(distribution(d), x)
quantile(d::PearsonAnalyticDistribution, p) = quantile(distribution(d), p)

"""
    distribution(d::PearsonAnalyticDistribution)

Returns a `Distributions.Distribution` object corresponding to `d`.
"""
function distribution end

"""
    relative_error(d::PearsonAnalyticDistribution, known)

Returns the (signed) relative error of `d`'s first non-fitted central moment and
the known coresponding central moment value `known`.  For `PearsonTypeIII` and
`PearsonTypeVI`, `known` should be the known fourth central moment.
"""
function relative_error end

# Thresholds

"""
    thresholds(d, nsigma::Real=3)

Compute upper and lower thresholds for distribution `d` that are equivalent to
`±nsigma` standard deviations of the standard normal distribution i.e.  `𝒩(µ=0,
σ=1)`.  `d` may be a `PearsonAnalyticDistribution`, a `PearsonTypeIV`, or a
`Distributions.Distribution`.  `nsigma` defaults to 3 if not given.
"""
function thresholds(pd::Union{PearsonAnalyticDistribution,PearsonTypeIV,Distribution}, nsigma::Real=3)
    nsigma = -abs(nsigma)
    p = cdf(Normal(), nsigma)
    quantile.(pd, (p, 1-p))
end
