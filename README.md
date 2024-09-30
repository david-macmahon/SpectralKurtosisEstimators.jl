# SpectralKurtosisEstimators.jl

This package provides the `SpectralKurtosisEstimator` type and methods for
working with them.  It is based on the series of spectral kurtosis papers by
G. M. Nita and D. E. Gary circa 2007-2016.  The `SpectralKurtosisEstimator` type
encapsulates the first four statstical moments of a spectral kurtosis estimator
using the generalized spectral kurtosis definition in equation 8 of:

> Nita, G. M. & Gary, D. E. [2010b] MNRAS 406, L60,
doi:10.1111/j.1745-3933.2010.00882.x.

The three defining parameters of a `SpectralKurtosisEstimator` are:

- `M`: number of pre-summed samples summed (e.g. "off-board")
- `N`: number of samples pre-summed (e.g "on-board")
- `d`: shape parameter for original voltage data
  - Use `1/2` for real single-pol inputs
  - Use `1` for complex single-pol inputs
  - Use `2` for complex dual-pol inputs(???)
    (maybe that should double `N` or `M` instead?)
