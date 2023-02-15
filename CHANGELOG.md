# CHANGELOG

## Pending

### Breaking changes

- [\#112](https://github.com/arkworks-rs/poly-commit/pull/112) Upgrade all dependencies to `0.4`.
- [\#82](https://github.com/arkworks-rs/poly-commit/pull/82) Argument `opening_challenge: F` for `open`,
  `check`, has been changed from `F` to `opening_challenges: &mut ChallengeGenerator`.

### Features

- [\#82](https://github.com/arkworks-rs/poly-commit/pull/82) Add multivariate opening challenge strategy. Integrate with sponge API.

### Improvements

### Bug fixes

## v0.3.0

### Breaking changes

- [\#78](https://github.com/arkworks-rs/poly-commit/pull/78) Fix `MarlinPC`'s `CommitterKey` to return the correct `supported_degree`.

### Features

### Improvements

### Bug fixes

## v0.2.0

- Initial release of `ark-poly-commit`.
