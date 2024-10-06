# SpectralKurtosisEstimators.jl

This package provides the `SKEstimator` type and methods for working with them.
It is based on the series of spectral kurtosis papers by G. M. Nita and D. E.
Gary circa 2007-2016:

* Radio Frequency Interference Excision Using Spectral‐Domain Statistics
  [Nita [2007]](https://doi.org/10.1086/520938)
* Statistics of the Spectral Kurtosis Estimator
  [Nita [2010a]](https://doi.org/10.1086/652409)
* A Wideband Spectrometer with RFI Detection
  [Gary [2010]](https://doi.org/10.1086/652410)
* The generalized spectral kurtosis estimator
  [Nita [2010b]](https://doi.org/10.1111/j.1745-3933.2010.00882.x)
* EOVSA Implementation of a Spectral Kurtosis Correlator for Transient Detection
  and Classification
  [Nita [2016]](http://doi.org/10.1142/S2251171716410099)

The `SKEstimator` type of this package encapsulates the first four statistical
moments of a *generalized spectral kurtosis estimator* using the definition from
equation 8 of [Nita [2010b]](https://doi.org/10.1111/j.1745-3933.2010.00882.x).
Other implementation details are from [Nita [2016]](
http://doi.org/10.1142/S2251171716410099).

## Overview

Spectral kurtosis is often used to identify outlier samples, especially in
dynamic power spectra.  This can be done by computing the spectral kurtosis
estimates from the input data and identifying values that lie outside some
lower/upper thresholds.  Computing the spectral kurtosis estimates involves only
relatively trivial arithmetic (summing, squaring, and dividing) so it is fairly
straightforward.  The trickier part is knowing what values to use for the
thresholds.

This package provides tools to model and analyze the statistical distributions
that arise from various parameterizations of the generalized spectral kurtosis
estimator.  It also provides a function, `skhat`, to compute the spectral
kurtosis estimates of a data array for a given set of spectral kurtosis
estimator parameters.

## SKEstimator

The generalized spectral kurtosis estimator is computed using nested sums.  The
inner sum has `N` addends while the outer sum has `M` addends.  The original
non-generalized estimator is a special case of the generalized estimator with
`N=1`.

The spectral kurtosis estimator also has a *shape* parameter, conventionally
referred to as `d`. This shape parameter is the same as the shape parameter of
the [gamma distribution](https://en.wikipedia.org/wiki/Gamma_distribution) of
the input samples (i.e. the addends of inner sum).  The value of `d` is half the
number of squared voltages that were added together per input sample.  Some
common values for `d` are shown in the table below, where `S` is the number of
squared voltages that were summed per input sample.

| Polarization | Voltages |  d  |
|:-------------|:---------|:---:|
| Single       | Real     | 1/2 |
| Stokes I     | Real     |  1  |
| Single       | Complex  |  1  |
| Stokes I     | Complex  |  2  |

Not surprisingly, `M`, `N`, and `d` are the three parameters of the
`SKEstimator` constructor:

    ske = SKEstimator(M, N, d)

where:

- `M`: number of outer sum addends
- `N`: number of inner sum addends
- `d`: half the number of squared voltages (pre-summed) per input sample

## Pearson distributions

[Pearson](https://en.wikipedia.org/wiki/Karl_Pearson) derived a [set of
distributions](https://en.wikipedia.org/wiki/Pearson_distribution) that
approximate an arbitrary distribution given only its first four (or three)
moments, which can be either computed from observations or derived analytically.
The papers listed above analytically derive formulas for the first four moments
of the spectral kurtosis estimator given the spectral kurtosis estimators `M`,
`N`, and `d`.  These moments can be used to create Pearson distributions of
various types.  This package supports Pearson Type IV, Pearson Type VI, and
Pearson Type III distributions.

### Pearson Type VI and Pearson Type III distributions

The Pearson Type VI distribution is actually a location shifted [beta prime
distribution]( https://en.wikipedia.org/wiki/Beta_prime_distribution).  The
Pearson Type III distribution is actually a location shifted [gamma
distribution]( https://en.wikipedia.org/wiki/Gamma_distribution).  In fact,
these now well known distributions originated from these earlier Pearson
distributions.  These Pearson distributions are represented internally by the
parameters of these well known distributions plus a location offset.

Because the Pearson Type VI and Pearson Type III distributions match only the
first three moments, their fourth moments will have an error relative to the
fourth moment of the distributions they are approximating.  Pearson
distributions with a lower magnitude relative error in the fourth moment better
fit the tails of the distribution they are approximating thereby making them
preferable for outlier detection.

Pearson Type VI and Pearson Type III distributions are represented by the
`PearsonTypeVI` and `PearsonTypeIII` types, resp.  They share a common abstract
super-type, `PearsonaAnalyticDistribution`, because they can be represented by
will known distributions supported by `Distributions.jl`.  The following
methods are supported for instances of these types:

- `mean(d::PearsonAnaliticType)` returns the mean of distribution `d`
- `var(d::PearsonAnaliticType)` returns the variance of distribution `d`
- `skewness(d::PearsonAnaliticType)` returns the skewness of distribution `d`
- `kurtosis(d::PearsonAnaliticType)` returns the excess kurtosis of distribution
  `d`
- `pdf(d::PearsonAnaliticType, x)` returns the probability density function of
  distribution `d` evaluated at `x`
- `cdf(d::PearsonAnaliticType, x)` returns the cumulative distribution function
of distribution `d` evaluated at `x`
- `quantile(d::PearsonAnaliticType, p)` returns `x` such that `cdf(d, x) == p`
  (i.e. the inverse cumulative distribution function)
- `thresholds(d::PearsonAnaliticType, nsigma)` returns the lower and upper
  values where the CDF of `d` equals the CDF of the standard normal
  distribution at `±nsigma`.
- `distribution(d::PearsonAnalyticType)` returns a `Distributions.jl`
  distribution corresponding to `d`.
- `relative_error(d::PearsonAnalyticType, u4)` returns the relative error
  between the fourth moment of `d` and `u4` (typically the fourth moment of an
  `SKEstimator`)
- `relative_error(d::PearsonAnalyticType, ske::SKEstimator)` returns
  `relative_error(d, ske.u4)`

### Pearson Type IV distribution

Pearson Type IV distributions match all of the first four moments.  It is only
possible to construct Pearson Type IV distributions if certain conditions are
met.  Pearson provided a formula, known as the *Pearson criterion*, that can be
used to determine whether the Pearson Type IV distribution may be constructed
for a given set of moments.

A Pearson Type IV distribution is represented by the `PearsonTypeIV` type.  The
following methods are supported for `PearsonTypeIV` instances:

- `mean(d::PearsonTypeIV)` returns the mean of distribution `d`
- `var(d::PearsonTypeIV)` returns the variance of distribution `d`
- `skewness(d::PearsonTypeIV)` returns the skewness of distribution `d`
- `kurtosis(d::PearsonTypeIV)` returns the excess kurtosis of distribution `d`
- `pdf(d::PearsonTypeIV, x)` returns the probability density function of
  distribution `d` evaluated at `x`

Notably missing are `distribution`, `relative_error`, and (for now) CDF related
functions `cdf`, `quantile`, and `thresholds`.

### Pearson criterion

The `pearson_criterion` function provided by this package computes the Pearson
criterion corresponding to a given `SKEstimator` (see equation 11 of [Nita
[2016]](http://doi.org/10.1142/S2251171716410099)).  The value of the Pearson
criterion, `κ`, determines the most suitable type(s) of Pearson distribution(s)
for a given set of moments, such as those for a given `SKEstimator`, as shown
here:

| Pearson criterion | Pearson distribution(s) |
|:-----------------:|:------------------------|
|      `κ < 0` .    | Type I (not supported)  |
|    `0 < κ < 1`    | Type VI                 |
|      `1 < κ`      | Type VI, Type III       |

## Computing spectral kurtosis estimates

Spectral kurtosis estimates can by computed by calling the `skhat` function,
which has these methods:

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
`s1` and `s2` on each call.  The returned Array will have the same number of
dimensions as `A`, but dimensions in `dims` will be 1.

All of these methods allocate the output Array each call.  In-place versions of
these methods do not yet exist, but broadcast can be used to store the result in
a suitable existing array without additional allocations.

    sk = similar(s1)
    sk .= skhat.(s1, s2, Ref(ske))
