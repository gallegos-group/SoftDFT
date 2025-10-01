# SoftDFT

**SoftDFT** is a modern, Julia-based codebase for modeling complex soft-matter systems using Ising density functional theory (iDFT) and its extensions.  
It is designed for high performance, readability, and extensibility.

---

## Overview

**Ising Density Functional Theory (iDFT)** provides a self-consistent framework for addressing these challenges by combining the statistical mechanics of internal state transitions (Ising-like models) with polymer DFT.  
This approach enables a unified description of spatially inhomogeneous, chemically responsive polymers wherein monomer ionization states and chain conformations coevolve under external stimuli.

SoftDFT implements iDFT and related frameworks in an efficient and extensible way, offering a platform for exploring charge regulation, polymer conformational transitions, and interfacial phenomena.

---

## Key Features

- **Julia-based**: readable, modern, high-performance language.
- **Modular architecture**: flexible representation of free energy functionals and molecular topologies.  
- **Advanced solvers**: Anderson mixing, nonlinear continuation methods.  
- **Extensible**: designed to be a foundation for new models and theoretical developments.  
- **Transparent & open-source**: encourages reproducibility and community contributions.  

---

## Installation

SoftDFT is written in Julia (≥ 1.9). To install:

```bash
# Clone the repository
git clone https://github.com/gallegos-group/SoftDFT.git
cd SoftDFT

# How to use section (in progress)

---

## Getting Started

# Example input files and workflows are provided in the `examples/` directory.  
# To run a simple test:

# ```julia
# include("examples/run_example.jl")
# ```

For details, see the [documentation](docs/) (in progress).

---

## Citation

If you use **SoftDFT** in academic work, please cite:

> Gallegos Group, *SoftDFT: A Julia-based platform for modeling soft matter with Ising density functional theory* (2025).  
> [https://github.com/gallegos-group/SoftDFT](https://github.com/gallegos-group/SoftDFT)

A `CITATION.cff` file is included in the repository for automatic citation formatting.

---

## License

This project is licensed under the [MIT License](LICENSE).

---

## Contributing

SoftDFT is under active development. Contributions are welcome via issues and pull requests.  
Please note: certain advanced modules remain under private development and are not included in this release.
