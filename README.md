# SpectralKurtosisEstimators.jl

This package provides the `SKEstimator` type and methods for working with them.
It is based on the series of spectral kurtosis papers by G. M. Nita and D. E.
Gary circa 2007-2016.  The `SKEstimator` type encapsulates the first four
statistical moments of a generalized spectral kurtosis estimator using the
generalized spectral kurtosis definition in equation 8 of:

> Nita, G. M. & Gary, D. E. [2010] MNRAS 406, L60,
doi:10.1111/j.1745-3933.2010.00882.x.

along with some implementation details from:

> Nita, G., Hickish, J., MacMahon, D., & Gary, D. 2016, JAI, 5, 1641009

## Overview

The three defining parameters of an `SKEstimator` are:

- `M`: number of pre-summed samples summed (e.g. "off-board")
- `N`: number of samples pre-summed (e.g "on-board")
- `d`: shape parameter for original voltage data
  - Use `1/2` for real voltages
  - Use `1` for complex voltages

Not surprisingly, these are the three parameters of the `SKEstimator`
constructor.

    ske = SKEstimator(M, N, d)

## Pearson probability curves

[Pearson](https://en.wikipedia.org/wiki/Karl_Pearson) derived a set of
[probability curves](https://en.wikipedia.org/wiki/Pearson_distribution) that
approximate the PDF of a distribution given its first four (or three) moments
derived from observations.  The papers listed above derive formulas for the
exact values of the first four moments of the spectral kurtosis estimator.
These values can be used to derive Pearson probability curves of various types.
This package supports Pearson Type III and Pearson type VI probability curves.
Pearson Type IV probability curves are not currently supported, but Pearson Type
VI curves can be used instead for practical applications (e.g. RFI detection)
with minimal error being introduced.

The Type III and Type VI probability curves do not use the fourth moment in
their construction, so the difference between their fourth moment and the
spectral kurtosis estimator's fourth moment can be used as an error metric of
the tails of the distributions.

The `pearson_type_iii` and `pearson_type_vi` functions return a `Distribution`
representing the Pearson probability curve of the corresponding type (III or VI)
for the given `SKEstimator` as well as the difference between the
`Distribution`'s fourth moment and the `SKEstimator`'s fourth moment expressed
as a relative error.

    skd_vi, u4err_vi = pearson_type_vi(ske)
    skd_iii, u4err_iii = pearson_type_iii(ske)

The returned `Distribution` objects can be used with other statistical functions
from the Julia ecosystem.  For example:

- `pdf(skd, x)` evaluates the probability density function of `skd` at `x`

- `cdf(skd, x)` evaluates the cumulative distribution function of `skd` at `x`

- `quantile(skd, q)` evaluates the inverse cumulative distribution function of
  `skd` at `q`.

## Thresholds

The *skewness* of the `SKEstimator` decreases relatively slowly, so the
distribution is asymmetric.  To get lower and upper thresholds for the Pearson
probability curve's distribution corresponding to the `±Nσ` thresholds of a
standard normal distribution, use the `SpectralKurtosisEstimators.thresholds`
function (shown here for `±3σ`):

    lower_vi, upper_vi = SpectralKurtosisEstimators.thresholds(skd_vi, 3)
    lower_iii, upper_iii = SpectralKurtosisEstimators.thresholds(skd_iii, 3)

## Computing spectral kurtosis estimates

Spectral kurtosis estimates can computed by calling the `skhat` function, which
has these methods:

    skhat(s1, s2, M, N=1, d=1)
    skhat(s1, s2, ske::SKEstimator)
    skhat(A, ske::SKEstimator; dims=:)

Compute the spectral kurtosis estimate from:

- `s1`: sum of power
- `s2`: sum of squared power
- `M`, `ske.M`: number of samples summed (e.g. "off-board")
- `N`, `ske.N`: number of samples pre-summed (e.g "on-board")
- `d`, `ske.d`: shape parameter for original voltage data
  - Use `1/2` for real voltages
  - Use `1` for complex voltages

In lieu of `s1` and `s2`, you can pass Array `A` and keyword argument `dims` to
have `s1` and `s1` be computed automatically.  Currently this method allocates
`s1` and `s2` on each call.

All of these methods allocate the output Array each call.  In-place versions of
these methods do not yet exist, but broadcast can be used to store the result in
a suitable existing array without additional allocations.

    sk = similar(s1)
    sk .= skhat.(s1, s2, Ref(ske))
