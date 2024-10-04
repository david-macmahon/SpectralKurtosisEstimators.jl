"""
Base type for all PearsonDistributions
"""
abstract type PearsonDistribution end

# Allow broadcasting over PearsonDistribution objects
Base.broadcastable(d::PearsonDistribution) = Ref(d)

"""
Supertype for types that can be represented by a known `Distribution`.  For
example, the `PearsonTypeIII` distribution is a location shifted `Gamma`
distribution and the `PearsonTypeVI` distribution is a location shifted
`BetaPrime` distribution.
"""
abstract type PearsonAnalyticDistribution <: PearsonDistribution end

"""
PearsonTypeIII distribution.  Essentially a location shifted `Gamma`
distribution.
"""
struct PearsonTypeIII <: PearsonAnalyticDistribution
    shape::Float64
    scale::Float64
    location::Float64
end

"""
    PearsonTypeIII(u2, u3)

Construct a PeasonTypeIII distribution having second and third central moments
given by `u2` and `u3`.
"""
function PearsonTypeIII(u2, u3)
    # The parameter naming convention of the 2010c SK paper differs from the
    # parameter naming convention of Distributions.jl for the Gamma
    # distribution.  The paper uses `β` for the *shape* parameter and `α` for
    # the *scale* parameter whereas `Distributions.Gamma` uses `α` for the shape
    # parameter and `θ` for the scale parameter.  To avoid (more) confusion, we
    # will refer to the shape and scale parameters here as `shape` and `scale`.
    # The *location* parameter called `δ` in the paper will be called `location`
    # here.
    shape = 4u2^3 / u3^2
    scale = u3 / 2u2
    location = 1 - 2u2^2 / u3

    PearsonTypeIII(shape, scale, location)
end

"""
    distribution(d::PearsonTypeIII)

Returns a `Distributions.Distribution` object corresponding to `PearsonTypeIII`
distribution `d`.
"""
function distribution(d::PearsonTypeIII)
    Gamma(d.shape, d.scale) + d.location
end

"""
    relative_error(pdiii::PearsonTypeIII, u4known)

Returns the (signed) relative error of `d`'s fourth central moment and the known
fourth central moment `u4known`.
"""
function relative_error(d::PearsonTypeIII, u4known)
    # Get local variables for pdiii's shape and scale fields
    (; shape, scale) = d
    u4dist = 3shape * (shape+2) * scale^4
    u4dist/u4known - 1
end

"""
PearsonTypeVI distribution.  Essentially a location shifted `BetaPrime`
distribution.
"""
struct PearsonTypeVI <: PearsonAnalyticDistribution
    a::Float64
    B::Float64
    location::Float64
end

"""
    PearsonTypeVI(u2, u3)

Construct a PeasonTypeVI distribution having second and third central moments
given by `u2` and `u3`.
"""
function PearsonTypeVI(u2, u3)
    radical = sqrt(16u2^4 + 4u3^2*u2 + u3^2)

    a = (
        32u2^5 - 4u3*u2^3 + 8u3^2*u2^2 + u3^2*u2 - u3^3
        + (8u2^3 - u3*u2 + u3^2) * radical
    ) / u3^3

    B = 2u2 * (4u2^2 + radical) / u3^2 + 3

    location = (B-a-1)/(B-1)

    PearsonTypeVI(a, B, location)
end

"""
    distribution(d::PearsonTypeVI)

Returns a `Distributions.Distribution` object corresponding to `PearsonTypeVI`
distribution `d`.
"""
function distribution(d::PearsonTypeVI)
    BetaPrime(d.a, d.B) + d.location
end

"""
Returns the (signed) relative error of `d`'s fourth central moment and the known
fourth central moment `u4known`.
"""
function relative_error(d::PearsonTypeVI, u4known)
    # Get local variables for d.a and d.B
    (; a, B) = d

    u4dist = (
        (3a * (a+B-1))
        /
        ((B-4) * (B-3) * (B-2) * (B-1)^4)
        *
        ((B+5)*a^2 + (B-1)*(B+5)*a + 2*(B-1)^2)
    )

    u4dist/u4known - 1
end

# Add statistical/distribution methods

mean(d::PearsonAnalyticDistribution) = mean(distribution(d))
var(d::PearsonAnalyticDistribution) = var(distribution(d))
skewness(d::PearsonAnalyticDistribution) = skewness(distribution(d))
kurtosis(d::PearsonAnalyticDistribution) = kurtosis(distribution(d))

pdf(d::PearsonAnalyticDistribution, x) = pdf(distribution(d), x)
cdf(d::PearsonAnalyticDistribution, x) = cdf(distribution(d), x)
quantile(d::PearsonAnalyticDistribution, p) = quantile(distribution(d), p)

# Thresholds

"""
    thresholds(d, nsigma::Real=3)

Compute upper and lower thresholds for distribution `d` that are equivalent to
`±nsigma` standard deviations of the standard normal distribution i.e.  `𝒩(µ=0,
σ=1)`.  `d` may be a `PearsonDistribution` or a `Distributions.Distribution`.
`nsigma` defaults to 3 if not given.
"""
function thresholds(pd::Union{PearsonDistribution,Distribution}, nsigma::Real=3)
    nsigma = -abs(nsigma)
    p = cdf(Normal(), nsigma)
    quantile.(pd, (p, 1-p))
end
