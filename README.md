# Structgen

LAMMPS structure generator for [LAMMPS](https://lammps.sandia.gov), developed by Vasilii Maksimov at the Functional Glasses and Materials Modeling Laboratory (FGM²L) at the University of North Texas, under the supervision of Dr. Jincheng Du.

**Structgen** is a command-line utility for generating initial atomic configurations for use with LAMMPS simulations of disordered systems, especially oxide glasses. It supports CHARGE-style atom definitions and allows flexible control over atomic counts, charges, and overall density.

## Features

- Command-line interface with `argparse`-based subcommands
- Generates LAMMPS `.structure` files with atomic positions and box boundaries
- Supports charged species via `atom_style charge`
- Automatically estimates simulation box size based on input density
- Modular design: atom placement strategies are handled via handler classes
- Written for reproducible initial structure generation workflows

## Installation

### Install from PyPI

Regular installation using pip:

```bash
pip install lammps-structgen
```

### Install directly from GitHub (bleeding-edge)

For the latest development version:

```bash
pip install git+https://github.com/superde1fin/structgen.git
```

### Clone & install locally

For local experimentation or contributing:

```bash
git clone https://github.com/superde1fin/structgen.git
cd structgen

# Standard install (no dev extras)
pip install .

# OR editable install (auto-reload while editing)
pip install -e .
```

> **Tip**: add the `[dev]` extra to pull in testing, linting, and publishing tools:
>
> ```bash
> pip install .[dev]
> ```

*Requires Python ≥ 3.9 and `numpy`.*

## Usage

After installation, invoke the CLI with:

```bash
structgen [subcommand] [options]
```

### Available subcommands:

- `charge`: Generate an initial structure for charge atom_style

### Example usage:

```bash
structgen charge -a 28.085 -a 16.00 2 -1.2 -d 2.2 -f 100 -o si_o_glass.structure
```

This generates a structure with 100 Si atoms (mass 28.085) and 200 O atoms (mass 16.00, charge -1.2), with an overall density of 2.2 g/cm³. The structure will be written to `si_o_glass.structure`.

## Documentation

All subcommands support `-h` or `--help` flags for detailed usage:

```bash
structgen charge --help
```

## Design Overview

The project is modular:

- `cli.py`: Main entry point, defines CLI and parses arguments
- `handlers/`: Contains one handler class per placement strategy (`CHARGE`, etc.)
- `constants.py`: Shared constants (default count and charge)
- `__init__.py`: Exposes version metadata
- `pyproject.toml`: Build and dependency metadata

Each handler:
- Parses user-provided atomic specs
- Calculates simulation box size from density
- Randomly places atoms while avoiding overlap
- Returns header and atom sections for LAMMPS input

## Development

Install dev dependencies:

```bash
pip install .[dev]
```

Run tests (if defined):

```bash
pytest
```

## License

GNU General Public License v3.0 (GPLv3). See [LICENSE](LICENSE) for details.

## Author

Vasilii Maksimov  
University of North Texas  
✉️ VasiliiMaksimov@my.unt.edu

## Links

- [LAMMPS Official Site](https://lammps.sandia.gov)
- [GitHub Repository](https://github.com/superde1fin/structgen)

