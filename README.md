# G-Functions

A small thermodynamics visualization project for exploring Gibbs energy functions and phase-boundary behavior in a binary A-B system.

## Overview

This project builds a phase diagram and visualizes how Gibbs energy curves and chemical potentials change with temperature. The workflow is separated into three parts:

- `g_functions/phase_diagram.py` — generates the equilibrium phase diagram data
- `g_functions/g_functions.py` — evaluates Gibbs energy and chemical potential functions
- `g_functions/plotting.py` — renders the interactive plots

## Project structure

- `g_functions/__init__.py`
- `g_functions/__main__.py`
- `g_functions/g_functions.py`
- `g_functions/phase_diagram.py`
- `g_functions/plotting.py`
- `GFunctions.py` — compatibility wrapper for the original script name
- `PhaseDiagram.py` — compatibility wrapper for the original script name
- `tests/` — validation for the model and dataset generation

## Installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Run the interactive plot

```bash
python -m g_functions
```

## Run the tests

```bash
pytest -q
```

## Notes

The original script style has been reorganized into a cleaner package layout without changing the underlying physics or plotting behavior. The goal is to make the project easier to understand, extend, and validate.
