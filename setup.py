from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="satay-cli",
    version="0.1.0",
    author="Anna Sintsova",
    author_email="ansintsova@ethz.ch",
    description="A tool for analyzing SATAY transposon insertion sequencing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    install_requires=[
        "click>=8.0.0",
        "pandas>=1.3.0",
        "pathlib>=1.0.1",
    ],
    entry_points={
        "console_scripts": [
            "satay=satay.cli:main",
        ],
    },
    package_data={
        "satay": ["py.typed"],  # If you're using type hints
    },
    include_package_data=True,
    zip_safe=False,
)
