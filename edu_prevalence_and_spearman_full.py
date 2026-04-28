"""
Generate (1) AD prevalence tables + heatmaps for:
  - EducationLevel × Age group
  - EducationLevel × Ethnicity
and (2) EducationLevel vs continuous covariates Spearman correlations:
  - CSV of rho/p
  - signed 1-row heatmap (blue=negative, red=positive)

Inputs (expected in same folder as this script):
  - alzheimers_disease_data_clean.csv

Outputs (written to same folder):
  - edu_x_agegrp_AD_prevalence.csv
  - edu_x_agegrp_AD_prevalence_heatmap.png
  - edu_x_ethnicity_AD_prevalence.csv
  - edu_x_ethnicity_AD_prevalence_heatmap.png
  - edu_assoc_spearman_educationlevel.csv
  - edu_assoc_spearman_signed_heatmap.png
"""

from __future__ import annotations

import os
import textwrap
from typing import List

import matplotlib

matplotlib.use("Agg")  # headless backend for reliable PNG export
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr


BASE_DIR = os.path.dirname(__file__)
DATA_CSV = os.path.join(BASE_DIR, "alzheimers_disease_data_clean.csv")

# Outputs
OUT_AGE_CSV = os.path.join(BASE_DIR, "edu_x_agegrp_AD_prevalence.csv")
OUT_AGE_PNG = os.path.join(BASE_DIR, "edu_x_agegrp_AD_prevalence_heatmap.png")
OUT_ETH_CSV = os.path.join(BASE_DIR, "edu_x_ethnicity_AD_prevalence.csv")
OUT_ETH_PNG = os.path.join(BASE_DIR, "edu_x_ethnicity_AD_prevalence_heatmap.png")

OUT_SPEAR_CSV = os.path.join(BASE_DIR, "edu_assoc_spearman_educationlevel.csv")
OUT_SPEAR_PNG = os.path.join(BASE_DIR, "edu_assoc_spearman_signed_heatmap.png")


def _wrap(s: str, width: int = 16) -> str:
    return textwrap.fill(str(s), width=width, break_long_words=False)


def load_data() -> pd.DataFrame:
    df = pd.read_csv(DATA_CSV)
    df.columns = [c.strip().lower() for c in df.columns]
    return df


def add_labels(df: pd.DataFrame) -> pd.DataFrame:
    # Age groups (<65, 65–74, >=75) to match your Stata do-file
    bins = [-np.inf, 65, 75, np.inf]
    age_labels = ["<65", "65–74", ">=75"]
    df["age_grp"] = pd.cut(df["age"], bins=bins, right=False, labels=age_labels)

    # Ethnicity mapping (from the notebook)
    eth_map = {0: "Caucasian", 1: "African American", 2: "Asian", 3: "Other"}
    if "ethnicity" in df.columns:
        df["ethnicity_lbl"] = df["ethnicity"].map(eth_map)

    # Education mapping (from the notebook)
    edu_map = {0: "None", 1: "High School", 2: "Bachelor's", 3: "Higher"}
    df["education_lbl"] = df["educationlevel"].map(edu_map)

    return df


def prevalence_table(
    df: pd.DataFrame,
    row: str,
    col: str,
    dx: str = "diagnosis",
    row_order: List[str] | None = None,
    col_order: List[str] | None = None,
) -> pd.DataFrame:
    sub = df[[row, col, dx]].dropna()
    prev = (
        sub.groupby([row, col], observed=False)[dx]
        .mean()
        .unstack(col)
        .astype(float)
        * 100.0
    )
    if row_order is not None:
        prev = prev.reindex(index=row_order)
    if col_order is not None:
        prev = prev.reindex(columns=col_order)
    return prev


def save_prevalence_heatmap(prev: pd.DataFrame, title: str, out_png: str) -> None:
    plt.close("all")
    fig_w = max(7.2, 1.2 * prev.shape[1] + 2.0)
    fig_h = max(4.8, 0.6 * prev.shape[0] + 2.4)
    fig, ax = plt.subplots(1, 1, figsize=(fig_w, fig_h), dpi=200)

    sns.heatmap(
        prev,
        ax=ax,
        cmap="Blues",  # light -> dark
        vmin=0,
        vmax=float(np.nanmax(prev.to_numpy())),
        annot=True,
        fmt=".1f",
        linewidths=0.5,
        linecolor="white",
        cbar=True,
        cbar_kws={"label": "AD prevalence (%)", "shrink": 0.9, "pad": 0.02},
    )
    ax.set_title(title, pad=12)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticklabels([_wrap(t.get_text(), 14) for t in ax.get_xticklabels()], rotation=0, ha="center")
    ax.set_yticklabels([t.get_text() for t in ax.get_yticklabels()], rotation=0)

    plt.tight_layout(rect=(0, 0, 0.98, 1))
    plt.savefig(out_png, dpi=200, bbox_inches="tight", pad_inches=0.25)


def spearman_educationlevel(df: pd.DataFrame) -> pd.DataFrame:
    cont_vars = [
        "age",
        "bmi",
        "alcoholconsumption",
        "physicalactivity",
        "dietquality",
        "sleepquality",
        "systolicbp",
        "diastolicbp",
        "cholesteroltotal",
        "cholesterolldl",
        "cholesterolhdl",
        "cholesteroltriglycerides",
        "mmse",
        "functionalassessment",
        "adl",
    ]

    rows = []
    for x in cont_vars:
        if x not in df.columns:
            continue
        s = df[["educationlevel", x]].dropna()
        rho, p = spearmanr(s["educationlevel"], s[x])
        rows.append(
            {
                "y": "educationlevel",
                "x": x,
                "N": int(len(s)),
                "rho": float(rho),
                "p": float(p),
                "abs_rho": float(abs(rho)),
            }
        )

    out = pd.DataFrame(rows).sort_values("rho")
    return out


def save_signed_spearman_heatmap(sp: pd.DataFrame, out_png: str) -> None:
    # 1-row heatmap, signed rho
    sp = sp.dropna(subset=["rho"]).sort_values("rho")
    series = sp.set_index("x")["rho"]
    cols = [_wrap(c, 16) for c in series.index.tolist()]
    data = pd.DataFrame([series.to_numpy()], index=[""], columns=cols)

    vmax = float(series.abs().max())
    fig_w = max(12, 0.55 * len(cols))
    fig_h = 3.2

    plt.close("all")
    plt.figure(figsize=(fig_w, fig_h), dpi=200)
    ax = sns.heatmap(
        data,
        cmap="RdBu_r",
        center=0.0,
        vmin=-vmax,
        vmax=vmax,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"shrink": 0.9, "pad": 0.03},
    )
    ax.set_title("EducationLevel vs continuous covariates (Spearman ρ)", pad=12)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_yticks([])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

    cbar = ax.collections[0].colorbar
    cbar.set_label("Spearman ρ\n(blue=neg, red=pos)", rotation=90, labelpad=18, fontsize=9)

    plt.tight_layout(rect=(0, 0, 0.86, 1))
    plt.savefig(out_png, dpi=200, bbox_inches="tight", pad_inches=0.6)


def main() -> None:
    df = add_labels(load_data())

    edu_order = ["None", "High School", "Bachelor's", "Higher"]
    age_order = ["<65", "65–74", ">=75"]
    eth_order = ["African American", "Asian", "Caucasian", "Other"]

    # 1) Prevalence: Education × Age group
    prev_age = prevalence_table(
        df,
        row="education_lbl",
        col="age_grp",
        row_order=edu_order,
        col_order=age_order,
    )
    prev_age.to_csv(OUT_AGE_CSV, index=True)
    save_prevalence_heatmap(prev_age, "AD prevalence by EducationLevel × Age group", OUT_AGE_PNG)

    # 2) Prevalence: Education × Ethnicity
    if "ethnicity_lbl" not in df.columns:
        raise RuntimeError("Missing ethnicity column in dataset; cannot build Ethnicity prevalence table.")
    prev_eth = prevalence_table(
        df,
        row="education_lbl",
        col="ethnicity_lbl",
        row_order=edu_order,
        col_order=eth_order,
    )
    prev_eth.to_csv(OUT_ETH_CSV, index=True)
    save_prevalence_heatmap(prev_eth, "AD prevalence by EducationLevel × Ethnicity", OUT_ETH_PNG)

    # 3) Spearman: EducationLevel vs continuous vars
    sp = spearman_educationlevel(df)
    sp.to_csv(OUT_SPEAR_CSV, index=False)
    save_signed_spearman_heatmap(sp, OUT_SPEAR_PNG)

    print("Wrote:")
    print(" -", OUT_AGE_CSV)
    print(" -", OUT_AGE_PNG)
    print(" -", OUT_ETH_CSV)
    print(" -", OUT_ETH_PNG)
    print(" -", OUT_SPEAR_CSV)
    print(" -", OUT_SPEAR_PNG)


if __name__ == "__main__":
    main()

