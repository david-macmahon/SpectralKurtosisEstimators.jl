# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project aims to adhere to [Semantic Versioning](https://semver.org/)
(with the usual 0.x convention that the minor version is the breaking boundary).

## [Unreleased]

## [0.4.0] - 2026-09-25

### Added
- Support for non-integer `M*N*d`: the moment formulas are valid for any
  `M*N*d > 0`, and the `isinteger(M*N*d)` restriction was lifted
- Parametric `SKEstimator{T}`: if `d` is a `Rational`, the moments are stored
  as exact `Rational{BigInt}` values; otherwise `T` is `Float64`
- Exact mixed-type constructor support (e.g. `SKEstimator(3, 2, 1//3)`)
- CI test workflow running the test suite on Julia 1.10 and the latest stable
  release
- Documentation of the new capabilities in the README and Getting Started page

### Changed
- `SKEstimator` is now a parametric type: `typeof` returns
  `SKEstimator{Float64}` or `SKEstimator{Rational{BigInt}}` rather than plain
  `SKEstimator`
- The constructors were simplified to a single positional method accepting any
  `Real` parameters plus a keyword-argument form
- `pearson_distribution` now throws an `ErrorException` for Pearson criterion
  values of 0 or less (previously it logged an error and returned `nothing`)

### Removed
- A leftover debug `try`/`catch` block from the `SKEstimator` constructor

## [0.3.3] - 2026-07-27

### Added
- Documenter.jl documentation
- Keyword-argument constructor for `SKEstimator`
- Doc string for the `SKEstimator` constructor

### Changed
- Relaxed the `Statistics` compat entry to 1.10
- Clarified documentation for `skhat(A, ske)`

### Fixed
- Typos in doc strings and README.md

## [0.3.2] - 2026-07-02

### Changed
- Bumped the `Roots` compat entry

## [0.3.1] - 2025-03-31

### Added
- `distribution(d::PearsonTypeIII, n)` returning the distribution of the
  average of `n` samples from a Pearson Type III distribution
- `std` method for Pearson distributions

### Fixed
- Typo in a doc string

## [0.3.0] - 2025-03-24

### Added
- `pearson_distribution` function selecting the most suitable Pearson
  distribution for an `SKEstimator`
- `cdf` and `quantile` support for `PearsonTypeIV`
- `[compat]` entries in Project.toml

### Removed
- Unused `loggammadiv` import

## [0.2.0] - 2025-03-16

### Changed
- Improved `skhat` for computing spectral kurtosis estimates of data Arrays
- Clarified README.md

## [0.1.0] - 2024-10-05

### Added
- `PearsonTypeIV` type
- `PearsonDistribution` type hierarchy (abstract supertypes for the Pearson
  distribution types)
- Split the Pearson Type III and Type VI implementations into separate files

## [0.0.2] - 2024-10-01

### Fixed
- Spelling of "criterion" and tweaks to doc strings and README.md

## [0.0.1] - 2024-10-01

### Added
- Initial release: the generalized spectral kurtosis estimator (`SKEstimator`)
  with the first four central moments from equation 9 of Nita & Gary (2010b)
- `mean`, `var`, `skewness`, and `kurtosis` methods for `SKEstimator`
- Pearson Type III and Type VI distributions with `mean`, `std`, `var`,
  `skewness`, `kurtosis`, `pdf`, `cdf`, `quantile`, and `thresholds` methods
- Pearson criterion computation with out-of-range warnings
- Unit tests and README

[Unreleased]: https://github.com/david-macmahon/SpectralKurtosisEstimators.jl/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/david-macmahon/SpectralKurtosisEstimators.jl/compare/v0.3.3...v0.4.0
[0.3.3]: https://github.com/david-macmahon/SpectralKurtosisEstimators.jl/compare/v0.3.2...v0.3.3
[0.3.2]: https://github.com/david-macmahon/SpectralKurtosisEstimators.jl/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/david-macmahon/SpectralKurtosisEstimators.jl/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/david-macmahon/SpectralKurtosisEstimators.jl/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/david-macmahon/SpectralKurtosisEstimators.jl/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/david-macmahon/SpectralKurtosisEstimators.jl/compare/v0.0.2...v0.1.0
[0.0.2]: https://github.com/david-macmahon/SpectralKurtosisEstimators.jl/compare/v0.0.1...v0.0.2
