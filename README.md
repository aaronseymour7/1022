# HFR_PACKAGE

A Python package for generating and analyzing **homodesmotic family reaction enthalpies** using RDKit, PuLP, and AaronTools.

---

## Features

- Balance Isogyric, Isodesmic, Hypohomodesmotic, and Homodesmotic reactions automatically.
- Compute reaction enthalpies from Gaussian output.
- References ATcT values for enthalpies of formation.
- Support for single and multiple computation workflows.
- Command-line interface (CLI) tools for ease of use.

---

## Installation

### Prerequisites

This package depends on some scientific libraries:

- [RDKit](https://www.rdkit.org/docs/Install.html)  
- [PuLP](https://coin-or.github.io/pulp/)  
- [AaronTools](https://github.com/QChASM/AaronTools.py) 

## Quick Start

### For Single Reaction
- `hfr [view,write,count] [Isogyric, Isodesmic, Hypohomodesmotic, Homodesmotic] [input SMILES or InChI] --m [method] --b [basis] --s [g,p,o]`
  Example:
    - `hfr write Isogyric 'InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3' --m B3LYP --b '6-31G(d)' --s g`

         writes an Isogyric reaction for butane with Gaussian opt+freq input files on B3LYP/6-31G(d) level of theory
  
  Files will be written in cwd with format `[R,P][INDEX]_[COEFF].[in]` e.g. (`R1_1.com`,`P2_3.in`)

To run jobs (from within reaction directory):
- `folder-run`

Once jobs have finished:
- `folder-compute`

Results will be in `opt_summary.txt` and summarized in `enthalpies_summary.csv`

### For Multiple Reactions

Create a line seperated plain Text file of format:

- `[1,2,3,4] [InChI]`
    1-Isogyric
    2-Isodesmic
    3-Hypohomodesmotic
    4-Homodesmotic
    space separation followed by InChi
Example:
```
1 'InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3'
3 'InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3'
```
To run reactions you must first

- ```multi-prep [your input file] --m [method] --b [basis] --s [g,p,o] ```



  
