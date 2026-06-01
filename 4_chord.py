#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Usage: python 4_chord.py [--input reports/category_map.csv] [--svg]

import argparse
import itertools
import os
import sys

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from matplotlib.colors import to_rgb


CATEGORIES = [
    "Aliphatic",
    "Ether",
    "Amide_Amine",
    "Alicyclic_Heterocycle",
    "Aromatic_Heteroaromatic",
    "Triazole_Click",
    "Rigid_Alkyne_Alkene",
    "Functionalised",
]

CATEGORY_COLORS = {
    "Aliphatic": "#F7C96A",
    "Ether": "#66B5D9",
    "Amide_Amine": "#66D1B0",
    "Functionalised": "#F7EE99",
    "Alicyclic_Heterocycle": "#C5D6FF",
    "Aromatic_Heteroaromatic": "#E88F55",
    "Triazole_Click": "#E3A9CA",
    "Rigid_Alkyne_Alkene": "#777777",
}


def build_cooccurrence(df, categories):
    cooc = pd.DataFrame(0, index=categories, columns=categories)
    for _, row in df.iterrows():
        present = [c for c in categories if row.get(c, False)]
        for a, b in itertools.combinations_with_replacement(present, 2):
            cooc.loc[a, b] += 1
            if a != b:
                cooc.loc[b, a] += 1
    return cooc


def _bezier(p0, p1, p2, steps=50):
    t = np.linspace(0, 1, steps)
    x = (1 - t) ** 2 * p0[0] + 2 * (1 - t) * t * p1[0] + t ** 2 * p2[0]
    y = (1 - t) ** 2 * p0[1] + 2 * (1 - t) * t * p1[1] + t ** 2 * p2[1]
    return x, y


def draw_chord(cooc, categories, out_html, out_svg=None):
    n = len(categories)
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    positions = {cat: (np.cos(theta[i]), np.sin(theta[i])) for i, cat in enumerate(categories)}

    mat = cooc.values.astype(float)
    max_val = mat.max()

    fig = go.Figure()

    for i, ci in enumerate(categories):
        for j, cj in enumerate(categories):
            if j <= i:
                continue
            val = mat[i, j]
            if val == 0:
                continue
            w = 3 + 20 * (np.log1p(val) / np.log1p(max_val))
            x0, y0 = positions[ci]
            x2, y2 = positions[cj]
            bx, by = _bezier((x0, y0), (0.0, 0.0), (x2, y2), steps=60)
            col_i = np.array(to_rgb(CATEGORY_COLORS[ci]))
            col_j = np.array(to_rgb(CATEGORY_COLORS[cj]))
            col_mid = (col_i + col_j) / 2
            color_str = f"rgb({col_mid[0]*255:.0f},{col_mid[1]*255:.0f},{col_mid[2]*255:.0f})"
            fig.add_trace(go.Scatter(
                x=list(bx), y=list(by),
                mode="lines",
                line=dict(width=w, color=color_str),
                opacity=0.45,
                hoverinfo="skip",
                showlegend=False,
            ))

    for cat in categories:
        x, y = positions[cat]
        node_size = 10 + 0.01 * cooc.sum(axis=1)[cat]
        fig.add_trace(go.Scatter(
            x=[x], y=[y],
            mode="markers+text",
            marker=dict(size=node_size, color=CATEGORY_COLORS[cat], line=dict(color="black", width=0.5)),
            text=[cat],
            textposition="bottom center",
            hoverinfo="text",
            showlegend=False,
        ))

    fig.update_layout(
        width=900, height=900,
        plot_bgcolor="white", paper_bgcolor="white",
        xaxis=dict(visible=False), yaxis=dict(visible=False),
        margin=dict(l=20, r=20, t=20, b=20),
    )

    fig.write_html(out_html, include_plotlyjs="cdn")
    print(f"Saved chord diagram → {out_html}")

    if out_svg:
        try:
            fig.write_image(out_svg)
            print(f"Saved SVG → {out_svg}")
        except Exception as e:
            print(f"[WARN] SVG export failed (install kaleido): {e}")


def main():
    ap = argparse.ArgumentParser(
        description="Interactive chord diagram of linker chemotype co-occurrence."
    )
    ap.add_argument("--input", default="reports/category_map.csv",
                    help="Category map CSV from linkeranalysis.py (default: reports/category_map.csv)")
    ap.add_argument("--outdir", default="graphs",
                    help="Output directory (default: graphs/)")
    ap.add_argument("--svg", action="store_true",
                    help="Also export a static SVG (requires kaleido)")
    args = ap.parse_args()

    if not os.path.exists(args.input):
        sys.exit(f"Input file not found: {args.input}\n"
                 "Run linkeranalysis.py first to generate reports/category_map.csv.")

    df = pd.read_csv(args.input)
    cooc = build_cooccurrence(df, CATEGORIES)

    os.makedirs(args.outdir, exist_ok=True)
    out_html = os.path.join(args.outdir, "chord_cooccurrence.html")
    out_svg = os.path.join(args.outdir, "chord_cooccurrence.svg") if args.svg else None

    draw_chord(cooc, CATEGORIES, out_html, out_svg)


if __name__ == "__main__":
    main()
