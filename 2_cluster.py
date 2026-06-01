#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Usage: python 2_cluster.py classified_linkers.xlsx

import argparse
import os
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from sklearn.decomposition import PCA


BASE_COLOR_MAP = {
    "Aliphatic": "#F7C96A",
    "Ether": "#66B5D9",
    "Amide_Amine": "#66D1B0",
    "Functionalised": "#F7EE99",
    "Alicyclic_Heterocycle": "#C5D6FF",
    "Aromatic_Heteroaromatic": "#E88F55",
    "Triazole_Click": "#E3A9CA",
    "Rigid_Alkyne_Alkene": "#777777",
    "Other": "#CCCCCC",
    "Hybrid": "#F5F5F5",
}


def morgan_fps(mols, radius=2, nbits=2048):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)
    return [gen.GetFingerprint(m) for m in mols]


def fps_to_numpy(fps):
    nbits = len(fps[0])
    X = np.zeros((len(fps), nbits), dtype=np.uint8)
    for i, fp in enumerate(fps):
        arr = np.zeros(nbits, dtype=np.int8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        X[i] = arr
    return X


def category_label(cat):
    if isinstance(cat, str) and cat.startswith("Hybrid(") and cat.endswith(")"):
        return "Hybrid"
    return cat if isinstance(cat, str) else "Other"


def main():
    ap = argparse.ArgumentParser(
        description="PCA of linker Morgan fingerprints → CSV for Prism."
    )
    ap.add_argument("xlsx", nargs="?", default="classified_linkers.xlsx",
                    help="Input Excel from 1_linkeranalysis.py (default: classified_linkers.xlsx)")
    ap.add_argument("--sheet", default="linkers")
    ap.add_argument("--radius", type=int, default=2)
    ap.add_argument("--nbits", type=int, default=2048)
    ap.add_argument("--output", default=None,
                    help="Output CSV path (default: <input_stem>_PCA.csv)")
    args = ap.parse_args()

    if not os.path.exists(args.xlsx):
        sys.exit(f"File not found: {args.xlsx}")

    df = pd.read_excel(args.xlsx, sheet_name=args.sheet)
    if "smiles" not in df.columns:
        sys.exit("Input file must contain a 'smiles' column.")

    mols, valid_idx = [], []
    for i, smi in enumerate(df["smiles"]):
        mol = Chem.MolFromSmiles(str(smi).strip())
        if mol is None:
            print(f"[WARN] Invalid SMILES at row {i}", file=sys.stderr)
            continue
        mols.append(mol)
        valid_idx.append(i)

    df_valid = df.iloc[valid_idx].reset_index(drop=True)

    fps = morgan_fps(mols, args.radius, args.nbits)
    X = fps_to_numpy(fps)
    pca = PCA(n_components=2, random_state=42)
    Xp = pca.fit_transform(X)

    var = pca.explained_variance_ratio_ * 100
    print(f"PC1: {var[0]:.1f}%  PC2: {var[1]:.1f}%  total: {var.sum():.1f}%")

    df_out = pd.DataFrame({
        "name": df_valid["name"],
        "smiles": df_valid["smiles"],
        "category": df_valid["category"],
        "cat_label": [category_label(c) for c in df_valid["category"]],
        "color": [BASE_COLOR_MAP.get(category_label(c), "#CCCCCC") for c in df_valid["category"]],
        "PC1": Xp[:, 0],
        "PC2": Xp[:, 1],
    })

    out_path = args.output or (os.path.splitext(args.xlsx)[0] + "_PCA.csv")
    df_out.to_csv(out_path, index=False)
    print(f"Saved PCA coordinates → {out_path}  ({len(df_out)} linkers)")


if __name__ == "__main__":
    main()
