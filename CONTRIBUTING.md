# Contributing to Bene Gesserit

Thank you for your interest in contributing to Bene Gesserit! This document provides guidelines and information for contributors.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Making Changes](#making-changes)
- [Testing](#testing)
- [Documentation](#documentation)
- [Pull Request Process](#pull-request-process)
- [Release Process](#release-process)

## Code of Conduct

This project adheres to a code of conduct. By participating, you are expected to uphold this code. Please report unacceptable behavior to the project maintainers.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally
3. Create a new branch for your feature or bugfix
4. Make your changes
5. Push to your fork and submit a pull request

## Development Setup

### Prerequisites

- Rust 1.70.0 or later
- Git
- A C compiler (for some dependencies)

### Setup Steps

```bash
# Clone your fork
git clone https://github.com/YOUR_USERNAME/bene-gesserit.git
cd bene-gesserit

# Install development dependencies
cargo install cargo-audit cargo-deny cargo-tarpaulin cargo-outdated

# Build the project
cargo build

# Run tests
cargo test

# Run clippy for linting
cargo clippy --all-targets --all-features

# Format code
cargo fmt
```

## Making Changes

### Branching Strategy

- `main` - Stable release branch
- `develop` - Development branch for upcoming releases
- `feature/*` - Feature branches
- `bugfix/*` - Bug fix branches
- `hotfix/*` - Critical fixes for production

### Commit Messages

Use conventional commit format:

```
type(scope): description

[optional body]

[optional footer]
```

Types: `feat`, `fix`, `docs`, `style`, `refactor`, `test`, `chore`

Examples:
- `feat(quantum): add ENAQT implementation`
- `fix(membrane): correct oscillation endpoint calculation`
- `docs(api): update BiologicalQuantumState documentation`

### Code Style

- Follow Rust naming conventions
- Use `rustfmt` for formatting (configuration in `rustfmt.toml`)
- Follow clippy suggestions (configuration in `clippy.toml`)
- Maximum line length: 100 characters
- Use meaningful variable and function names
- Add documentation comments for public APIs

### Scientific Code Guidelines

- Ensure all equations are properly documented with references
- Include units in variable names or comments where applicable
- Validate physical constraints and biological realism
- Add comprehensive tests for numerical accuracy
- Document mathematical derivations in comments or docs

## Testing

### Test Types

1. **Unit Tests**: Test individual functions and modules
2. **Integration Tests**: Test component interactions
3. **Doctests**: Test code examples in documentation
4. **Benchmarks**: Performance testing for critical paths

### Running Tests

```bash
# Run all tests
cargo test

# Run tests with coverage
cargo tarpaulin --all-features

# Run benchmarks
cargo bench

# Run specific test
cargo test test_biological_quantum_state

# Run tests with output
cargo test -- --nocapture
```

### Test Guidelines

- Write tests for all public APIs
- Include edge cases and error conditions
- Use property-based testing for mathematical functions
- Test numerical accuracy with known solutions
- Mock external dependencies

## Documentation

### Documentation Types

1. **API Documentation**: Rustdoc comments on public items
2. **Examples**: Working code examples
3. **Tutorials**: Step-by-step guides
4. **Reference**: Comprehensive API reference

### Writing Documentation

- Use rustdoc syntax for code documentation
- Include examples in doc comments
- Document mathematical concepts and biological context
- Link to relevant papers and references
- Keep documentation up-to-date with code changes

### Building Documentation

```bash
# Build documentation
cargo doc --all-features --no-deps

# Build and open documentation
cargo doc --all-features --no-deps --open

# Test documentation examples
cargo test --doc
```

## Pull Request Process

### Before Submitting

1. Ensure all tests pass
2. Run clippy and fix warnings
3. Format code with rustfmt
4. Update documentation if needed
5. Add or update tests for your changes
6. Update CHANGELOG.md

### PR Requirements

- Clear description of changes
- Link to relevant issues
- Include test results
- Follow commit message conventions
- Ensure CI passes
- Get review from maintainer

### Review Process

1. Automated checks (CI, tests, linting)
2. Code review by maintainers
3. Address feedback and update PR
4. Final approval and merge

## Release Process

### Version Numbering

We follow [Semantic Versioning](https://semver.org/):
- MAJOR: Incompatible API changes
- MINOR: Backward-compatible functionality additions
- PATCH: Backward-compatible bug fixes

### Release Steps

1. Update version in `Cargo.toml`
2. Update CHANGELOG.md
3. Create release branch
4. Test thoroughly
5. Create GitHub release
6. Publish to crates.io (automated)

## Getting Help

- Check existing issues and documentation
- Ask questions in GitHub discussions
- Join our community chat (if available)
- Contact maintainers directly for sensitive issues

## Scientific Contributions

### Biological Accuracy

- Ensure biological parameters are physiologically realistic
- Cite relevant literature for biological processes
- Validate against experimental data where possible
- Consider cellular constraints and energy budgets

### Mathematical Rigor

- Provide mathematical derivations for new algorithms
- Validate numerical stability
- Test conservation laws where applicable
- Document approximations and their limitations

### Performance Considerations

- Profile critical code paths
- Optimize for typical biological simulation scales
- Consider memory usage for large systems
- Document computational complexity

Thank you for contributing to Bene Gesserit! 