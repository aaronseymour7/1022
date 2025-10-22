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
`hfr [view,write,count] [Isogyric, Isodesmic, Hypohomodesmotic, Homodesmotic] [input SMILES or InChI] --m [method] --b [basis] --s [g,p,o]`

|      | Software   |
|:----:|:-----------|
| g    | Gaussian   |
| p    | Psi4       |
| o    | ORCA       |

  **Example:**
     `hfr write Isogyric 'InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3' --m B3LYP --b '6-31G(d)' --s g`

  - writes an Isogyric reaction for butane with Gaussian opt+freq input files on B3LYP/6-31G(d) level of theory
  
  **NOTE** Files will be written in cwd with format `[R,P][INDEX]_[COEFF].[in]` e.g. (`R1_1.com`,`P2_3.in`)

To run jobs (from within reaction directory):
 `folder-run`

Once jobs have finished:
 `folder-compute`

Results will be in `opt_summary.txt` and summarized in `enthalpies_summary.csv`

### For Multiple Reactions

Create a **line-separated plain text** file of the format:
 
 `[1,2,3,4] [InChI]`
    
 Where the first int term indicates the **reaction type**:

|      | Reaction Type          |
|:----:|:-----------------------|
| 1    | Isogyric               |
| 2    | Isodesmic              |
| 3    | Hypohomodesmotic       |
| 4    | Homodesmotic           |

Each line should contain a **space-separated** pair

**Example**:
```
1 'InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3'
3 'InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3'
```
To run reactions you must first

 `multi-prep [your input file] --m [method] --b [basis] --s [g,p,o] `

  This will create a directory [int].mhfr for each input molecule. Inside each directory will be the reactant and product input files along with an index.txt that stores reactant and product names. 

To run all reactions you must then run:

`multi-run`

This will submit all unique jobs.

Once all jobs are complete, run:

`multi-compute`

All result data can be found in `enthalpies_summaries.csv` as well as `opt_summary.txt`. For individual reaction data, go to the `.mhfr` directory corresponding to the reaction of interest's line number in the input file.

## Optional Enhancements

### Generate Conformers

When using the input.txt multi-reaction input format, after prepping: 

`multi-prep [your input file] --m [method] --b [basis] --s [g,p,o] `

run

`multi-mc --rmsd [rmsd]`

- This will generate and run conformers where the quantity is a function of rotatable bonds and rmsd is the overlap similarity metric. default RMSD = .15

After jobs are finished, run: 

`multi-ec`

`multi-compute`

### Single Point Energy Calculation Add-on

After running the standard workflow Single Point Energy calculations can be run using:

`multi-sprun --m [method] --b [basis] --s [g,p,o]`

|      | Software   |
|:----:|:-----------|
| g    | Gaussian   |
| p    | Psi4       |
| o    | ORCA       |

This will create a `unique_files/spec/` directory filled with Single Point input files for the desired level of theory in the desired software and submit the jobs.

Once the jobs are finished, run

`multi-spcomp`

This will compute enthalpies of formation using the Single Point energy values. Additionally, each `.mhfr` file will have a subdirectory `/spec/` added with its reaction specific SP energy values. Results can be found in `sp_enthalpies_summary.csv`



  
