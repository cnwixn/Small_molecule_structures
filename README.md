# Small molecule structures

Turn a SMILES string into a 2D depiction (SVG) and a 3D conformer (SDF).

## What it does

- **2D** — draws the molecule with RDKit CoordGen and saves a clean black-and-white SVG (`<name>_2D.svg`). Bonds are 2 px, atom labels use the BW palette, so the output drops straight into a figure.
- **3D** — adds hydrogens, embeds a conformer with ETKDGv3, optimises it with the MMFF94 force field and writes an SDF (`<name>_3D.sdf`).
- Prints the molecular formula, atom count and bond count.

## Install

```
pip install -r requirements.txt
```

## Usage

```
python smiles_to_structure.py
python smiles_to_structure.py "CCO"
python smiles_to_structure.py --outdir results --name aspirin "CC(=O)OC1=CC=CC=C1C(=O)O"
```

Options:

| Option | Default | Meaning |
| --- | --- | --- |
| `smiles` | acetylsalicylic acid | SMILES string of the molecule |
| `-o`, `--outdir` | `.` | output directory |
| `-n`, `--name` | `molecule` | base name of the output files |
| `--seed` | `42` | random seed for the 3D embedding |
| `--width`, `--height` | `500`, `300` | SVG canvas size in pixels |

## Sketching the 2D structure

Draw the molecule in [Ketcher](https://github.com/epam/ketcher) (ketcher-standalone 3.7.0 or newer), export the SMILES string, and pass it to the script.

## Notes

- The 3D embedding is deterministic. Two runs with the same seed give identical coordinates, so figures and structures are reproducible.
- To look at the 3D structure, open the SDF in PyMOL, ChimeraX or any other molecular viewer. A notebook is no longer needed.
