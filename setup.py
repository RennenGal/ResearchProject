"""
Setup configuration for the Protein Data Collector package.
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="protein-data-collector",
    version="0.1.0",
    author="Protein Data Collector Team",
    author_email="team@protein-data-collector.com",
    description="A bioinformatics application for collecting TIM barrel protein data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/your-org/protein-data-collector",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=7.4.0",
            "pytest-asyncio>=0.21.0",
            "hypothesis>=6.80.0",
            "pytest-cov>=4.1.0",
            "black>=23.0.0",
            "isort>=5.12.0",
            "flake8>=6.0.0",
            "mypy>=1.4.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "protein-collector=protein_data_collector.cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "protein_data_collector": ["config/*.json", "sql/*.sql"],
    },
)