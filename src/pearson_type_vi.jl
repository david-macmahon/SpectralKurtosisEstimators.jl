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

# Documented in pearson_distributions.jl
function distribution(d::PearsonTypeVI)
    BetaPrime(d.a, d.B) + d.location
end

# Documented in pearson_distributions.jl
function relative_error(d::PearsonTypeVI, u4known::Real)
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
