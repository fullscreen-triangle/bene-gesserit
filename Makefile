# Makefile for Bene Gesserit Development

.PHONY: help build test check fmt clippy doc bench clean install-tools setup ci

# Default target
help:
	@echo "Bene Gesserit Development Commands"
	@echo "=================================="
	@echo ""
	@echo "Setup:"
	@echo "  setup          - Install all development tools"
	@echo "  install-tools  - Install cargo tools"
	@echo ""
	@echo "Development:"
	@echo "  build          - Build the project"
	@echo "  test           - Run all tests"
	@echo "  check          - Run cargo check"
	@echo "  fmt            - Format code"
	@echo "  clippy         - Run clippy lints"
	@echo "  doc            - Build documentation"
	@echo "  bench          - Run benchmarks"
	@echo ""
	@echo "Quality:"
	@echo "  ci             - Run all CI checks locally"
	@echo "  audit          - Security audit"
	@echo "  deny           - Check with cargo-deny"
	@echo "  coverage       - Generate test coverage"
	@echo ""
	@echo "Maintenance:"
	@echo "  clean          - Clean build artifacts"
	@echo "  update         - Update dependencies"

# Setup development environment
setup: install-tools
	@echo "🔧 Setting up development environment..."
	cargo build
	@echo "✅ Setup complete!"

# Install development tools
install-tools:
	@echo "📦 Installing development tools..."
	cargo install cargo-audit cargo-deny cargo-tarpaulin cargo-outdated cargo-edit
	@echo "✅ Tools installed!"

# Build the project
build:
	@echo "🔨 Building project..."
	cargo build --all-features

# Run tests
test:
	@echo "🧪 Running tests..."
	cargo test --all-features

# Run cargo check
check:
	@echo "🔍 Checking code..."
	cargo check --all-features

# Format code
fmt:
	@echo "🎨 Formatting code..."
	cargo fmt --all

# Run clippy
clippy:
	@echo "📎 Running clippy..."
	cargo clippy --all-targets --all-features -- -D warnings

# Build documentation
doc:
	@echo "📚 Building documentation..."
	cargo doc --all-features --no-deps --open

# Run benchmarks
bench:
	@echo "⚡ Running benchmarks..."
	cargo bench --all-features

# Run all CI checks locally
ci: fmt clippy test doc
	@echo "🚀 Running CI checks..."
	cargo audit
	cargo deny check
	@echo "✅ All CI checks passed!"

# Security audit
audit:
	@echo "🔒 Running security audit..."
	cargo audit

# Check with cargo-deny
deny:
	@echo "🚫 Running cargo-deny..."
	cargo deny check

# Generate test coverage
coverage:
	@echo "📊 Generating test coverage..."
	cargo tarpaulin --all-features --out Html

# Clean build artifacts
clean:
	@echo "🧹 Cleaning build artifacts..."
	cargo clean

# Update dependencies
update:
	@echo "🔄 Updating dependencies..."
	cargo update
	cargo outdated

# Quick development cycle
dev: fmt clippy test
	@echo "✅ Development cycle complete!"

# Release preparation
release-prep: ci coverage
	@echo "📦 Preparing for release..."
	@echo "✅ Release preparation complete!" 