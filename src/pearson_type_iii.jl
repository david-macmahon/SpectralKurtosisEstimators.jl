"""
PearsonTypeIII distribution.  Essentially a location shifted `Gamma`
distribution.
"""
struct PearsonTypeIII <: PearsonAnalyticDistribution
    shape::Float64
    scale::Float64
    location::Float64
    n::Int
end

"""
    PearsonTypeIII(u2, u3)

Construct a PeasonTypeIII distribution having second and third central moments
given by `u2` and `u3`.
"""
function PearsonTypeIII(u2, u3, n=1)
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

    PearsonTypeIII(shape, scale, location, n)
end

# Documented in pearson_distributions.jl
function distribution(d::PearsonTypeIII, n::Integer=d.n)
    Gamma(n * d.shape, d.scale)/n + d.location
end

# Documented in pearson_distributions.jl
function relative_error(d::PearsonTypeIII, u4known::Real)
    # TODO Incorporate d.n into this calculation
    # Get local variables for d.shape and d.scale
    (; shape, scale) = d
    u4dist = 3shape * (shape+2) * scale^4
    u4dist/u4known - 1
end
