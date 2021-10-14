# Changelog

This is the Changelog of MathUtils.
This file should be kept up to date following [these guidelines](https://keepachangelog.com/en/1.0.0/)

## [Unreleased]

### Added


### Changed


### Fixed

## [1.8.1] - 2021-10-12
### Fixed
- Bug in Interp1d Eval on std::vector

## [1.8] - 2021-10-11

### Major changes
- Interp1d derived class added to reintroduce extrapolation and saturation behavior 
    - Interp1dLinearSaturate
    - Interp1dLinearExtrapolate

## [1.7.2] - 2021-09-15

### Added
- Bessel functions derivatives implemented


## [1.7.1] - 2021-09-14

### Added
- First kind Hankel function implemented and tested

## [1.7] - 2021-08-12

### Major changes
- Permissive option removed from Interp1d and LookUpTable

### Fixed
- Boost URL changed

## [1.5] - 2021-04-01
`><(((°>      ><(((°>       ><(((°>       ><(((°>`

### Major changes
- Boost::boost removed from MathUtils dependencies
- MathUtils components requiring boost, in distinct MathUtilsBoost.h include file

### Added
- NEW CHANGELOG !!!
- googletest for unit testing, not applied on tests, only one example in test_BoostFunctions.cpp