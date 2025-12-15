#!/usr/bin/env python3
"""
Setup script for Essential Gene Analysis package.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README
readme_path = Path(__file__).parent / "README.md"
long_description = readme_path.read_text() if readme_path.exists() else ""

setup(
    name="essential-gene-analysis",
    version="1.0.0",
    author="Roghaye Mohammadbeygi",
    author_email="rbm5893@psu.edu",
    description="Tool for identifying conditionally essential genes and rescue reactions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/essential-gene-analysis",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "cobra>=0.26.0",
        "pandas>=1.3.0",
        "numpy>=1.20.0",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "pyyaml>=5.4.0",
        "tqdm>=4.60.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0.0",
            "pytest-cov>=2.12.0",
            "black>=21.0.0",
            "flake8>=3.9.0",
            "mypy>=0.900",
        ],
        "vis": [
            "matplotlib-venn>=0.11.0",
            "upsetplot>=0.6.0",
            "networkx>=2.6.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "essential-genes=src.cli:main",
        ],
    },
    include_package_data=True,
    package_data={
        "": ["*.yaml", "*.yml", "*.json"],
    },
)
