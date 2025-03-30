"""
PearsonTypeIV distribution.  The CDF of this function is computed using a
function that approximates the real CDF to (nearly) machine precision over a
finite but practical range.
"""
struct PearsonTypeIV{m,ν,a,λ} <: PearsonApproximatedDistribution where {m,ν,a,λ}
    u2::Float64
    u3::Float64
    u4::Float64
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

    PearsonTypeIV{m,ν,a,λ}(u2, u3, u4)
end

# Add statistical/distribution methods

mean(d::PearsonTypeIV) = 1.0
std(d::PearsonTypeIV) = sqrt(var(d))
var(d::PearsonTypeIV) = d.u2
skewness(d::PearsonTypeIV) = d.u3 / d.u2^(3/2)
kurtosis(d::PearsonTypeIV) = d.u4 / d.u2^2 - 3

function _pearsontypeiv_pdf(m, ν, a, λ, x)
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

"""
    pdf(d::PearsonTypeIV, x)

Return the PDF of `d` evaluated at `x`.
"""
function pdf(::PearsonTypeIV{m,ν,a,λ}, x) where {m,ν,a,λ}
    _pearsontypeiv_pdf(m, ν, a, λ, x)
end

"""
    _cdflimits(d::PearsonTypeIV, z=eps()) -> (lo, hi)

Return lower and upper limits of the CDF of `d` suitable for use with
`_cdffun` (as well as `cdf` and `quantile`).  These limits are determined by
solving for the values below and above the mean of `d` where the PDF equals
`z`.
"""
@memoize function _cdflimits(d::PearsonTypeIV, z=eps())
    µ = mean(d)
    lo = find_zero(x->pdf(d, x)-z, (-Inf, µ))
    hi = find_zero(x->pdf(d, x)-z, (µ, Inf))
    lo, hi
end

"""
    _cdffun(d::PearsonTypeIV) -> ApproxFun.Fun
    _cdffun(d::PearsonTypeIV, lo, hi) -> ApproxFun.Fun

Return approximated function `f(x)` that will return the approximate (to nearly
machine precision) CDF of `d` at `x`.  If `lo` and `hi` are given, they set the
domain of the approximated CDF function of `d`.  If they are not given, suitable
values are determined automatically.  The returned function is only valid from
`lo` to `hi`.  For values outside that range the returned function will return
`0.0`.  NB: Values above `hi` will return a CDF value of `0.0` rather than
`1.0`!
"""
@memoize function _cdffun(d::PearsonTypeIV, lo, hi)
    cumsum(Fun(x->pdf(d, x), lo..hi))
end

function _cdffun(d::PearsonTypeIV)
    lo, hi = _cdflimits(d)
    _cdffun(d, lo, hi)
end

"""
    cdf(d::PearsonTypeIV, x)
    cdf(d::PearsonTypeIV, x, lo, hi)

Return the approximated (to nearly machine precision) CDF of `d` evaluated at
`x`.  If `lo` and `hi` are given, they set the domain of the approximated CDF
function of `d`.  If they are not given, suitable values are determined
automatically.  If `x` is less than `lo`, `0.0` is returned.  If `x` is greater
than `hi`, `1.0` is returned.
"""
function cdf(d::PearsonTypeIV, x, lo, hi)
    if x < lo
        return 0.0
    elseif x > hi
        return 1.0
    else
        return _cdffun(d, lo, hi)(x)
    end
end

function cdf(d::PearsonTypeIV, x)
    lo, hi = _cdflimits(d)
    cdf(d, x, lo, hi)
end

"""
    quantile(d::PearsonTypeIV, p) -> x
    quantile(d::PearsonTypeIV, p, lo, hi) -> x

Returns the inverse CDF of `d` evaluated at `p` by solving this equation for
`x`:

    cdf(d, x) == p

The CDF is computed using an approximated function (see `cdf(d::PearsonTypeIV,
x, lo, hi)`).  If `lo` and `hi` are given, they set the domain of the
approximated CDF function.  If they are not given, suitable values are
determined automatically.  The return value for different values of `p` is shown
here:

| `p`                             | Return value                       |
|:-------------------------------:|:-----------------------------------|
| `p < 0`                         | N/A (throws error)                 |
| `0 <= p < cdf(d, lo)`           | `lo`                               |
| `cdf(d, lo) <= p <= cdf(d, hi)` | `x` such that `cdf(d, x) == p`     |
| `cdf(d, hi) < p <= 1`           | `hi`                               |
| `1 < p`                         | N/A (throws error)                 |
"""
function quantile(d::PearsonTypeIV, p, lo, hi)
    cdffun = _cdffun(d, lo, hi)
    if p < 0.0 || p > 1.0
        @error "p must be in interval [0,1]"
    elseif p < cdffun(lo)
        lo
    elseif p <= cdffun(hi)
        roots(cdffun - p) |> only
    else
        hi
    end
end

function quantile(d::PearsonTypeIV, p)
    lo, hi = _cdflimits(d)
    quantile(d::PearsonTypeIV, p, lo, hi)
end
