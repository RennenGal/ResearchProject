.PHONY: install test clean

install:
	pip install -e .
	pip install -r requirements-dev.txt

test:
	pytest

clean:
	rm -rf build/ dist/ *.egg-info/ .pytest_cache/ .mypy_cache/ htmlcov/ coverage.xml .coverage
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete
