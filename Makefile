.PHONY: help install install-dev test test-unit test-property test-integration lint format type-check clean docs run-collector

# Default target
help:
	@echo "Available targets:"
	@echo "  install       - Install package and dependencies"
	@echo "  install-dev   - Install package with development dependencies"
	@echo "  test          - Run all tests"
	@echo "  test-unit     - Run unit tests only"
	@echo "  test-property - Run property-based tests only"
	@echo "  test-integration - Run integration tests only"
	@echo "  lint          - Run code linting"
	@echo "  format        - Format code with black and isort"
	@echo "  type-check    - Run type checking with mypy"
	@echo "  clean         - Clean build artifacts and cache"
	@echo "  docs          - Build documentation"
	@echo "  run-collector - Run the protein collector CLI"

# Installation targets
install:
	pip install -e .

install-dev:
	pip install -e ".[dev]"
	pip install -r requirements-dev.txt
	pre-commit install

# Testing targets
test:
	pytest

test-unit:
	pytest -m "unit"

test-property:
	pytest -m "property"

test-integration:
	pytest -m "integration"

test-coverage:
	pytest --cov=protein_data_collector --cov-report=html --cov-report=term

# Code quality targets
lint:
	flake8 protein_data_collector/ tests/
	black --check protein_data_collector/ tests/
	isort --check-only protein_data_collector/ tests/

format:
	black protein_data_collector/ tests/
	isort protein_data_collector/ tests/

type-check:
	mypy protein_data_collector/

# Utility targets
clean:
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	rm -rf htmlcov/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

docs:
	@echo "Documentation build not yet implemented"

# Application targets
run-collector:
	python -m protein_data_collector.cli --help

# Database setup (placeholder)
setup-db:
	@echo "Database setup not yet implemented"

# Development server (placeholder)
dev-server:
	@echo "Development server not yet implemented"