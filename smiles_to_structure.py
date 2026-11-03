#!/usr/bin/env python3
"""Turn a SMILES string into a 2D depiction and a 3D conformer.

Usage:
    python smiles_to_structure.py
    python smiles_to_structure.py "CCO"
    python smiles_to_structure.py --outdir results --seed 7 "CC(=O)OC1=CC=CC=C1C(=O)O"

Without a SMILES argument the script falls back to the built-in example,
acetylsalicylic acid.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

DEFAULT_SMILES = "CC(=O)OC1=CC=CC=C1C(=O)O"
DEFAULT_SEED = 42


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert a SMILES string into a 2D SVG depiction and a 3D SDF conformer."
    )
    parser.add_argument(
        "smiles",
        nargs="?",
        default=DEFAULT_SMILES,
        help="SMILES string of the molecule (default: acetylsalicylic acid)",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        default=".",
        help="directory for the output files (default: current directory)",
    )
    parser.add_argument(
        "-n",
        "--name",
        default="molecule",
        help="base name for the output files (default: molecule)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help="random seed for the 3D embedding, keeps runs reproducible (default: 42)",
    )
    parser.add_argument(
        "--width", type=int, default=500, help="SVG width in pixels (default: 500)"
    )
    parser.add_argument(
        "--height", type=int, default=300, help="SVG height in pixels (default: 300)"
    )
    return parser.parse_args(argv)


def molecule_from_smiles(smiles: str) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        sys.exit(f"Could not parse SMILES: {smiles}")
    return mol


def write_2d_svg(mol: Chem.Mol, path: Path, width: int, height: int) -> Path:
    """Draw a clean black-and-white 2D depiction and save it as SVG."""
    mol_2d = Chem.RemoveHs(mol)
    mol_2d.RemoveAllConformers()
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol_2d)

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 2.0
    opts.minFontSize = 14
    opts.additionalAtomLabelPadding = 0.1
    opts.useBWAtomPalette()
    opts.prepareMolsBeforeDrawing = True

    drawer.DrawMolecule(mol_2d)
    drawer.FinishDrawing()
    path.write_text(drawer.GetDrawingText(), encoding="utf-8")
    return path


def write_3d_sdf(mol: Chem.Mol, path: Path, seed: int) -> Path:
    """Embed a 3D conformer with ETKDG, optimise it with MMFF94 and save it as SDF."""
    mol_3d = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol_3d, params) != 0:
        sys.exit("3D embedding failed for this molecule")
    AllChem.MMFFOptimizeMolecule(mol_3d)

    writer = Chem.SDWriter(str(path))
    try:
        writer.write(mol_3d)
    finally:
        writer.close()
    return path


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    mol = molecule_from_smiles(args.smiles)
    svg_path = write_2d_svg(mol, outdir / f"{args.name}_2D.svg", args.width, args.height)
    sdf_path = write_3d_sdf(mol, outdir / f"{args.name}_3D.sdf", args.seed)

    print(f"SMILES:      {args.smiles}")
    print(f"Formula:     {rdMolDescriptors.CalcMolFormula(mol)}")
    print(f"Atoms:       {mol.GetNumAtoms()}")
    print(f"Bonds:       {mol.GetNumBonds()}")
    print(f"2D SVG:      {svg_path}")
    print(f"3D SDF:      {sdf_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
