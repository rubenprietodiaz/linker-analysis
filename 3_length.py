#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Usage: python 3_length.py classified_linkers.xlsx

import argparse
import os
import sys
import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem


N_CONFS_DEFAULT = 30
MAX_ITERS = 200


def generate_conformers(smiles, n_confs):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        cids = AllChem.EmbedMultipleConfs(mol, n_confs, params=params)
        not_converged = 0
        for cid in cids:
            result = AllChem.UFFOptimizeMolecule(mol, confId=cid, maxIters=MAX_ITERS)
            if result == 1:
                not_converged += 1
        if not_converged:
            print(f"[WARN] {not_converged}/{len(cids)} conformers did not converge: {smiles[:60]}",
                  file=sys.stderr)
        return mol, list(cids)
    except Exception:
        return None


def max_heavy_atom_distance(conf, heavy_idx):
    max_d = 0.0
    n = len(heavy_idx)
    for i in range(n):
        pi = conf.GetAtomPosition(heavy_idx[i])
        for j in range(i + 1, n):
            d = pi.Distance(conf.GetAtomPosition(heavy_idx[j]))
            if d > max_d:
                max_d = d
    return max_d


def _process_row(args):
    smi, n_confs = args
    out = generate_conformers(smi, n_confs)
    if out is None:
        return np.nan
    mol, cids = out
    if not cids:
        return np.nan
    heavy_idx = [i for i in range(mol.GetNumAtoms()) if mol.GetAtomWithIdx(i).GetAtomicNum() > 1]
    return float(max(max_heavy_atom_distance(mol.GetConformer(cid), heavy_idx) for cid in cids))


def main():
    ap = argparse.ArgumentParser(
        description="Maximum end-to-end distance via 30-conformer ETKDGv3 + UFF ensemble."
    )
    ap.add_argument("xlsx", nargs="?", default="classified_linkers.xlsx",
                    help="Input Excel file from 1_linkeranalysis.py (default: classified_linkers.xlsx)")
    ap.add_argument("--sheet", default="linkers",
                    help="Sheet name (default: linkers)")
    ap.add_argument("--output", default=None,
                    help="Output Excel path (default: <input_stem>_distances.xlsx)")
    ap.add_argument("--nconfs", type=int, default=N_CONFS_DEFAULT,
                    help=f"Conformers per molecule (default: {N_CONFS_DEFAULT})")
    ap.add_argument("--workers", type=int, default=min(cpu_count(), 8),
                    help="Parallel worker processes (default: min(cpu_count, 8))")
    args = ap.parse_args()

    if not os.path.exists(args.xlsx):
        sys.exit(f"File not found: {args.xlsx}")

    df = pd.read_excel(args.xlsx, sheet_name=args.sheet)
    if "smiles" not in df.columns:
        sys.exit("Input file must contain a 'smiles' column.")

    out_path = args.output or (os.path.splitext(args.xlsx)[0] + "_distances.xlsx")
    tasks = [(str(smi), args.nconfs) for smi in df["smiles"]]

    results = []
    with Pool(processes=args.workers) as pool:
        for res in tqdm(pool.imap(_process_row, tasks), total=len(tasks), desc="Conformers"):
            results.append(res)

    df_out = df.copy()
    df_out["dist_max_A"] = results
    df_out.to_excel(out_path, index=False)
    print(f"Done. Output saved to: {out_path}")


if __name__ == "__main__":
    main()
