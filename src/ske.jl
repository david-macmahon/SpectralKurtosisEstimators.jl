struct SKEstimator{T<:Real}
    M::Int
    N::Int
    d::T

    u1::T
    u2::T
    u3::T
    u4::T
end

"""
    SKEstimator(M, N=1, d=1)
    SKEstimator(; M, N=1, d=1)

Construct a generalized spectral kurtosis estimator with `M` outer sum addends,
`N` inner sum addends, and shape parameter `d`.  The first four central moments
(`u1`, `u2`, `u3`, `u4`) of the estimator are precomputed using the formulas
from equation 9 of:

> "Monthly Notices of the Royal Astronomical Society". 406, L60-L64 (2010)
> doi:10.1111/j.1745-3933.2010.00882.x

Arguments:

- `M`: number of outer sum addends (must be `>= 2`)
- `N`: number of inner sum addends (must be `>= 1`, defaults to `1`)
- `d`: half the number of squared voltages summed together per input sample
  (must be `> 0`, defaults to `1`).  This is the same as the shape parameter of
  the [gamma distribution](https://en.wikipedia.org/wiki/Gamma_distribution).
  - Use `1/2` for single-pol real voltages
  - Use `1` for single-pol complex voltages or Stokes I from real voltages
  - Use `2` for Stokes I from complex voltages

`M` and `N` must be integers, but `d` may be any positive real number.  The
moment formulas are a closed-form evaluation of the gamma-function expressions
of equation 9 of the reference, valid for any `M*N*d > 0`.  The eltype `T` of
the returned `SKEstimator{T}` is `Rational{BigInt}` when `d` is a `Rational`
(in which case the moments are exact), and `Float64` otherwise.  The moments
can be accessed via `mean`, `var`, `skewness`, and `kurtosis`, or directly
through the `u1`, `u2`, `u3`, and `u4` fields.

See also: [`skhat`](@ref), [`pearson_distribution`](@ref).
"""
function SKEstimator(M::Real, N::Real=1, d::Real=1)
    # Rational d gets exact moments in Rational{BigInt} arithmetic; anything
    # else (including plain Integers, for backwards-compatible pragmatics)
    # computes in Float64
    T = d isa Rational ? Rational{BigInt} : Float64

    isinteger(M) || error("value of M ($M) must be an integer")
    isinteger(N) || error("value of N ($N) must be an integer")
    Mint, Nint = Int(M), Int(N)

    Mint >= 2 || error("M ($Mint) must be >= 2")
    Nint >= 1 || error("N ($Nint) must be >= 1")
    d > 0 || error("d ($d) must be > 0")

    # Do all arithmetic in T: Float64 avoids Int64 overflow in M^2 etc. for
    # huge M, and Rational{BigInt} keeps Rational inputs exact
    M, N, d = T(M), T(N), T(d)
    Nd = N*d
    MNd = M*Nd

    u1=one(T)

    # Central moments, equation 9 of Nita & Gary (2010).  The gamma-function
    # ratios of equation 9 are evaluated in closed (product) form, which is
    # algebraically identical for all real MNd > 0:
    #   Γ(MNd+2)/Γ(MNd+4) = 1/((MNd+2)(MNd+3))
    #   Γ(MNd+2)/Γ(MNd+6) = 1/((MNd+2)(MNd+3)(MNd+4)(MNd+5))
    #   Γ(MNd+2)/Γ(MNd+8) = 1/((MNd+2)(MNd+3)(MNd+4)(MNd+5)(MNd+6)(MNd+7))
    u2=(
        (2M^2 * Nd * (1+Nd))
        /
        ((M-1) * (6 + 5M*Nd + M^2*Nd^2))
    )

    u3=(
        (8M^3 * Nd * (Nd+1))
        * ((-2 + Nd * (-5 + M *(4 + Nd))))
        / ((M-1)^2 * (MNd+5) * (MNd+4) * (MNd+3) * (MNd+2))
    )

    u4=(
        (12M^4 * Nd * (Nd+1))
        * (24 + Nd*(48 + 84Nd + M*(-32 + Nd*(-245 - 93Nd + M*(125 + Nd*(68 + M + (3 + M)*Nd))))))
        /
        ((M-1)^3 * (MNd+7) * (MNd+6) * (MNd+5) * (MNd+4) * (MNd+3) * (MNd+2))
    )

    SKEstimator{T}(Mint, Nint, d, u1, u2, u3, u4)
end

SKEstimator(; M, N=1, d::Real=1) = SKEstimator(M, N, d)

# Statistics for SKEstimator
mean(ske::SKEstimator) = ske.u1
var(ske::SKEstimator) = ske.u2
skewness(ske::SKEstimator) = ske.u3 / sqrt(ske.u2)^3
kurtosis(ske::SKEstimator) = ske.u4 / ske.u2^2 - 3

"""
    pearson_criterion(ske) -> Real

Return the Pearson criterion for SKEstimator `ske`.
"""
function pearson_criterion(ske::SKEstimator)
    pearson_criterion(ske.u2, ske.u3, ske.u4)
end

"""
    pearson_distribution(ske) -> PearsonDistribution

Construct the optimal PearsonDistribution for SKEstimator `ske`.  If the Pearson
criterion for `ske` is between 0 and 1, a `PearsonTypeIV` distribution will be
returned.  If the Pearson criterion is 1 or greater, the `PearsonTypeIII` or
`PearsonTypeVI` distribution for `ske` with the lower relative error in the
fourth moment will be returned.  A Pearson criterion of 0 or less is not
supported and will throw an `ErrorException`.
"""
function pearson_distribution(ske::SKEstimator)
    κ = pearson_criterion(ske)
    # TODO Verify the correct thing to do when κ == 1
    if 0 < κ < 1
        PearsonTypeIV(ske)
    elseif 1 <= κ
        err3 = relative_error(PearsonTypeIII, ske)
        err6 = relative_error(PearsonTypeVI, ske)
        if err3 < err6
            PearsonTypeIII(ske)
        else
            PearsonTypeVI(ske)
        end
    else
        error("Pearson criterion $κ <= 0 is not supported")
    end
end

# Construct specific PearsonDistributions from an SKEstimator
PearsonTypeIII(ske::SKEstimator) = PearsonTypeIII(ske.u2, ske.u3)
PearsonTypeIV(ske::SKEstimator) = PearsonTypeIV(ske.u1, ske.u2, ske.u3, ske.u4)
PearsonTypeVI(ske::SKEstimator) = PearsonTypeVI(ske.u2, ske.u3)

# Relative error functions for PearsonDistribution `d` and SKEstimator `ske`
relative_error(d::PearsonAnalyticDistribution, ske::SKEstimator) = relative_error(d, ske.u4)

# Relative error for PearsonAnalyticDistribution type `D` for SKEstimator `ske`
relative_error(D::Type{<:PearsonAnalyticDistribution}, ske::SKEstimator) = relative_error(D(ske), ske.u4)
