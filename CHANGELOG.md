# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.3] - 2026-02-23

### Deprecated
- This crate is deprecated in favor of `complex-bessel`, a pure Rust implementation.

## [0.1.2] - 2024-12-19

### Added
- Scaled function variants: `J_scaled`, `Y_scaled`, `I_scaled`, `K_scaled`, `Ai_scaled`, `Bi_scaled`
- Comprehensive test suite for all functions
- Platform support documentation (macOS, Windows MSVC/GNU, Linux)

### Changed
- Renamed "Advanced API" to "Low-level API" for clarity
- Enhanced README with detailed scaling factor explanations
- Improved API documentation with zeta parameter clarification

### Deprecated
- None

### Removed
- None

### Fixed
- None

### Security
- None

## [0.1.1] - 2024-12-19

### Added
- MSVC compiler support for Windows builds
- CHANGELOG.md file

### Changed
- None

### Deprecated
- None

### Removed
- None

### Fixed
- Build compatibility issues on Windows with MSVC compiler

### Security
- None

## [0.1.0] - 2024-12-19

### Added
- Initial release with complex Bessel functions and Airy functions
- Simple API for single value calculations
- Advanced API for multiple value calculations
- Thread-safe implementation
- Automatic C binding generation 