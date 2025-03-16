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
the [gamma distribution](https://en.wikipedia.org/wiki/Gamma_distribution)
(i.e. half the number of squares summed into each sample).  In signal processing
terms, `d` is half the number of power values (i.e. squared voltage values) that
were summed together per input sample.  Some common values for `d` are shown in
the table here:

| Polarization | Voltages |  d  | Summation per sample                        |
|:-------------|:---------|:---:|:--------------------------------------------|
| Single       | Real     | 1/2 | `V^2`                                       |
| Stokes I     | Real     |  1  | `Vx^2 + Vy^2`                               |
| Single       | Complex  |  1  | `Re(V)^2 + Im(V)^2`                         |
| Stokes I     | Complex  |  2  | `Re(Vx)^2 + Im(Vx)^2 + Re(Vy)^2 + Im(Vy)^2` |

Not surprisingly, `M`, `N`, and `d` are the three parameters of the
`SKEstimator` constructor:

    ske = SKEstimator(M, N, d)

where:

- `M`: number of outer sum addends
- `N`: number of inner sum addends
- `d`: half the number of squared voltages summed together per input sample (aka
  the *shape* parameter).

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
well known distributions supported by `Distributions.jl`.  The following
methods are supported for instances of these types:

- `mean(d::PearsonAnalyticDistribution)` returns the mean of distribution `d`
- `var(d::PearsonAnaliyicDistribution)` returns the variance of distribution `d`
- `skewness(d::PearsonAnalyticDistribution)` returns the skewness of
  distribution `d`
- `kurtosis(d::PearsonAnalyticDistribution)` returns the excess kurtosis of
  distribution `d`
- `pdf(d::PearsonAnalyticDistribution, x)` returns the probability density
  function of distribution `d` evaluated at `x`
- `cdf(d::PearsonAnalyticDistribution, x)` returns the cumulative distribution
  function of distribution `d` evaluated at `x`
- `quantile(d::PearsonAnalyticDistribution, p)` returns `x` such that `cdf(d, x)
  == p` (i.e. the inverse cumulative distribution function)
- `thresholds(d::PearsonAnalyticDistribution, nsigma)` returns the lower and
  upper values where the CDF of `d` equals the CDF of the standard normal
  distribution at `±nsigma`.
- `distribution(d::PearsonAnalyticDistribution)` returns a `Distributions.jl`
  distribution corresponding to `d`.
- `relative_error(d::PearsonAnalyticDistribution, u4)` returns the relative
  error between the fourth moment of `d` and `u4` (typically the fourth moment
  of an `SKEstimator`)
- `relative_error(d::PearsonAnalyticDistribution, ske::SKEstimator)` returns
  `relative_error(d, ske.u4)`

### Pearson Type IV distribution

Pearson Type IV distributions match all of the first four moments.  It is only
possible to construct Pearson Type IV distributions if certain conditions are
met.  Pearson provided a formula, known as the *Pearson criterion* (see below),
that can be used to determine whether a Pearson Type IV distribution may be
constructed from a distribution's first four moments.  It is worth noting that
when `Nd` (i.e. the product of spectral kurtosis estimator parameters `N` and
`d`) is greater than 14 it is not possible to construct a Pearson Type IV
distribution.

A Pearson Type IV distribution is represented by the `PearsonTypeIV` type.  The
following methods are supported for `PearsonTypeIV` instances:

- `mean(d::PearsonTypeIV)` returns the mean of distribution `d`
- `var(d::PearsonTypeIV)` returns the variance of distribution `d`
- `skewness(d::PearsonTypeIV)` returns the skewness of distribution `d`
- `kurtosis(d::PearsonTypeIV)` returns the excess kurtosis of distribution `d`
- `pdf(d::PearsonTypeIV, x)` returns the probability density function of
  distribution `d` evaluated at `x`

Understandably missing from that list are `distribution` (`PearsonTypeIV` has no
corresponding distribution from `Distributions.jl`) and `relative_error`
(`PearsonTypeIV` has no fourth moment error, by definition, but in theory
`relative_error` could return the error in the fifth moment).  More glaringly
missing are CDF related functions `cdf`, `quantile`, and `thresholds`, which
will be added in a future version.

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
|    `0 < κ < 1`    | Type IV                 |
|      `1 < κ`      | Type VI, Type III       |

## Computing spectral kurtosis estimates

The `skhat` function provides various methods for computing spectral kurtosis
estimates.

### Spectral kurtosis estimates from pre-calculated `s1` and `s2`

Spectral kurtosis estimates can by computed from pre-calculated `s1` and `s2` by
calling one of these `skhat` methods:

    skhat(s1, s2, M, N=1, d=1)
    skhat(s1, s2, ske::SKEstimator)

where:

- `s1`: sum of power
- `s2`: sum of squared power
- `M`, `ske.M`: number of samples summed (e.g. "off-board" or "outer sum")
- `N`, `ske.N`: number of samples pre-summed (e.g "on-board" or "inner sum")
- `d`, `ske.d`: half the number of squared voltages (pre-summed) per input
  sample
  - Use `1/2` for single-pol real voltages
  - Use `1` for single-pol complex voltages or Stokes I from real voltages
  - Use `2` for Stokes I from complex voltages

These methods can be used with broadcast to store the output into a suitably
sized pre-allocated Array.

### Spectral kurtosis estimates of data in an Array

The spectral kurtosis estimates of an Array cam be computed by calling this
`skhat` method:

    skhat(A, ske::SKEstimator; dims=ndims(A))

This computes the generalized spectral kurtosis estimate of `A` along `dims` as
specified by the `M` and `N` fields of `ske`.  `dims` must be an integer and
defaults to the last dimension of `A`.  `ske.d` should be set to half the number
of real samples that were pre-summed into each element of `A`.

Named tuple `(; s1, sk)` is returned, where `s1` is the summed power and `sk` is
the spectral kurtosis.

`s1` and `sk` will have the same dimensions as `A` except that dimension `dims`
will be `size(A, dims) ÷ (M*N)`.  If `M*N` does not divide `size(A, dims)`
evenly then some number (less than `M*N`) of samples from the end of the `dims`
dimension of `A` will not be used.

This method allocate the intermediate and output Arrays each call.
