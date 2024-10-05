"""
PearsonTypeIV distribution.  The CDF of this function is computed using a
function that approximates the real CDF to (nearly) machine precision over a
finite but practical range.
"""
struct PearsonTypeIV <: PearsonApproximatedDistribution
    u2::Float64
    u3::Float64
    u4::Float64

    r::Float64
    m::Float64
    ν::Float64
    a::Float64
    λ::Float64
end

"""
    PearsonTypeIV(u2, u3, u4)

Construct a PeasonTypeIV distribution having second, third, and fourth central
moments given by `u2`, `u3`, and `u4`.
"""
function PearsonTypeIV(u1, u2, u3, u4)
    pc = pearson_criterion(u2,u3,u4)
    if !(0 < pc < 1)
        error("Pearson criterion $pc is not in (0..1)")
    end

    B1 = u3^2 / u2^3
    B2 = u4 / u2^2
    r = 6*(B2-B1-1) / (2B2-3B1-6)

    m = (r+2) / 2
    ν = -r*(r-2)*sqrt(B1) / sqrt(16*(r-1) - B1*(r-2)^2)
    a = sqrt(u2*(16*(r-1) - B1*(r-2)^2)) / 4
    λ = u1 - (r-2)*sqrt(u2*B1)/4

    PearsonTypeIV(u2, u3, u4, r, m, ν, a, λ)
end

# Add statistical/distribution methods

mean(d::PearsonTypeIV) = 1.0
var(d::PearsonTypeIV) = d.u2
skewness(d::PearsonTypeIV) = d.u3 / d.u2^(3/2)
kurtosis(d::PearsonTypeIV) = d.u4 / d.u2^2 - 3

function pdf(d::PearsonTypeIV, x)
    # Get local variables for some of d's fields
    (; m, ν, a, λ) = d

    #= # Non-log formula
    real(
        (gamma(complex(m, ν/2)) * gamma(complex(m, -ν/2)))
        /
        (a * sqrt(π) * gamma(m-1/2) * gamma(m))
        *
        ((1+((x - λ)/a)^2)^-m * exp(-ν * atan(x-λ, a)))
    )
    =#

    real(exp(
        -log(a) - log(π)/2 + loggamma(complex(m, ν/2)) + loggamma(complex(m, -ν/2))
        - loggamma(m-1/2) - loggamma(m) - m * log1p(((x-λ)/a)^2) - ν * atan(x-λ, a)
    ))
end
