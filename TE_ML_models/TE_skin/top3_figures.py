#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import re

from matplotlib.colors import TwoSlopeNorm

# Set up args
parser = argparse.ArgumentParser()
parser.add_argument("--model", choices=["control", "benzonase"], required=True,
    help="Which source treatment to train on (and which output folder/model name to use).")
args = parser.parse_args()

# Set variables
MODEL = args.model
OUTPUT = f"{MODEL}_output"

print("Plotting Top-3 inclusion heatmaps...")

def donor_key(label: str):
    s = str(label).strip()
    m = re.search(r'(\d+)$', s)
    return (0, int(m.group(1))) if m else (1, s)

def topk_inclusion_matrix_fixed_axes(df, ylabels, xlabels, topk_col="Top3", sep=";"):
    row_index = {lab: i for i, lab in enumerate(ylabels)}
    col_index = {lab: j for j, lab in enumerate(xlabels)}

    M = np.zeros((len(ylabels), len(xlabels)), dtype=float)
    counts = df["True_Label"].astype(str).value_counts().to_dict()

    for _, r in df.iterrows():
        t = str(r["True_Label"]).strip()
        if t not in row_index:
            continue
        topk = [x.strip() for x in str(r[topk_col]).split(sep) if x.strip()]
        i = row_index[t]
        denom = counts.get(t, 0)
        if denom == 0:
            continue
        for pred in topk:
            if pred in col_index:
                j = col_index[pred]
                M[i, j] += 1.0 / denom   # "fraction of samples where pred appears in Top-3"

    return M

def collect_label_union(df_list, topk_col="Top3", sep=";"):
    yset = set()
    xset = set()
    for df in df_list:
        yset |= set(df["True_Label"].astype(str).str.strip().unique().tolist())
        # include everything that appears in Top3 (and also true labels)
        xset |= set(df["True_Label"].astype(str).str.strip().unique().tolist())
        for s in df[topk_col].astype(str):
            for item in s.split(sep):
                item = item.strip()
                if item:
                    xset.add(item)
    ylabels = sorted(yset, key=donor_key)
    xlabels = sorted(xset, key=donor_key)

    return ylabels, xlabels

def plot_heatmap(M, ylabels, xlabels, title, outpath, vmax):
    plt.figure(figsize=(10, 8))
    im = plt.imshow(M, aspect="auto", interpolation="nearest", cmap="Blues", vmin=0, vmax=vmax)
    ax = plt.gca()
    ax.set_xticks(np.arange(-0.5, len(xlabels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(ylabels), 1), minor=True)
    ax.grid(which="minor", linestyle="-", linewidth=0.5)
    plt.colorbar(im, label="Fraction of samples where donor appears in Top-3")
    plt.yticks(np.arange(len(ylabels)), ylabels)
    plt.xticks(np.arange(len(xlabels)), xlabels, rotation=90)
    plt.ylabel("True label")
    plt.xlabel("Top-3 predicted label (inclusion)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, bbox_inches="tight")
    plt.close()

# Load both
df_control = pd.read_csv(f"{OUTPUT}/source_{MODEL}_sink_control_predictions.csv")
df_treat   = pd.read_csv(f"{OUTPUT}/source_{MODEL}_sink_treatment_predictions.csv")

# Shared axes (union of labels)
ylabels, xlabels = collect_label_union([df_control, df_treat], topk_col="Top3", sep=";")

# Build matrices on the same axes
M_control = topk_inclusion_matrix_fixed_axes(df_control, ylabels, xlabels, topk_col="Top3", sep=";")
M_treat   = topk_inclusion_matrix_fixed_axes(df_treat,   ylabels, xlabels, topk_col="Top3", sep=";")

# Shared vmax
vmax = float(max(M_control.max(), M_treat.max()))

# Plot both using the same vmax
plot_heatmap(
    M_control, ylabels, xlabels,
    title=f"Top-3 inclusion heatmap (source: {MODEL}, sink: control)",
    outpath=f"{OUTPUT}/{MODEL}_top3_inclusion_control.pdf",
    vmax=vmax
)

plot_heatmap(
    M_treat, ylabels, xlabels,
    title=f"Top-3 inclusion heatmap (source: {MODEL}, sink: benzonase)",
    outpath=f"{OUTPUT}/{MODEL}_top3_inclusion_treatment.pdf",
    vmax=vmax
)

print("Wrote:",
      f"{OUTPUT}/{MODEL}_top3_inclusion_control.pdf",
      "and",
      f"{OUTPUT}/{MODEL}_top3_inclusion_treatment.pdf")

D = M_treat - M_control
v = float(np.max(np.abs(D)))

norm = TwoSlopeNorm(vmin=-v, vcenter=0.0, vmax=v)

plt.figure(figsize=(10, 8))
im = plt.imshow(D, aspect="auto", interpolation="nearest", cmap="PuOr", norm=norm)
ax = plt.gca()
ax.set_xticks(np.arange(-0.5, len(xlabels), 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(ylabels), 1), minor=True)
ax.grid(which="minor", linestyle="-", linewidth=0.5)
plt.colorbar(im, label="Δ inclusion (treatment − control)")
plt.yticks(np.arange(len(ylabels)), ylabels)
plt.xticks(np.arange(len(xlabels)), xlabels, rotation=90)
plt.ylabel("True label")
plt.xlabel("Top-3 predicted label (inclusion)")
plt.title(f"Δ Top-3 inclusion heatmap (treatment − control) (source: {MODEL})")
plt.tight_layout()
plt.savefig(f"{OUTPUT}/{MODEL}_top3_inclusion_delta.pdf", bbox_inches="tight")
plt.close()

# Per-donor recall at 3 plot
print("Plotting Top-3 recall barplot...")

def topk_hit(true_label, topk_str, sep=";"):
    items = [x.strip() for x in str(topk_str).split(sep) if x.strip()]
    return int(str(true_label).strip() in items)

def recall_at_k_by_class(df, topk_col="Top3", sep=";"):
    df = df.copy()
    df["hit"] = [topk_hit(t, k, sep=sep) for t, k in zip(df["True_Label"], df[topk_col])]
    out = df.groupby("True_Label")["hit"].mean()
    return out

ctrl = pd.read_csv(f"{OUTPUT}/source_{MODEL}_sink_control_predictions.csv")
treat = pd.read_csv(f"{OUTPUT}/source_{MODEL}_sink_treatment_predictions.csv")

r3_ctrl = recall_at_k_by_class(ctrl, "Top3")
r3_treat = recall_at_k_by_class(treat, "Top3")

labels = sorted(set(r3_ctrl.index) | set(r3_treat.index), key=donor_key)
x = np.arange(len(labels))
w = 0.4

plt.figure(figsize=(12, 4))
plt.bar(x - w/2, [r3_ctrl.get(l, np.nan) for l in labels], width=w, label="Control")
plt.bar(x + w/2, [r3_treat.get(l, np.nan) for l in labels], width=w, label="Treatment")
plt.xticks(x, labels, rotation=90)
plt.ylabel("Top-3 hit rate")
plt.title(f"Per-donor Top-3 hit rate (source: {MODEL})")
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUTPUT}/{MODEL}_recall_at3_per_donor_control_vs_treat.pdf", bbox_inches="tight")
plt.close()

# top 1 vs 3 vs 5
print("Plotting Top-1 vs Top-3 vs Top-5 barplot...")

def topk_accuracy(df, col="Top3"):
    return np.mean([topk_hit(t, k) for t, k in zip(df["True_Label"], df[col])])

top1_ctrl = np.mean(ctrl["True_Label"].astype(str).str.strip() == ctrl["Predicted_Label"].astype(str).str.strip())
top1_treat = np.mean(treat["True_Label"].astype(str).str.strip() == treat["Predicted_Label"].astype(str).str.strip())
top3_ctrl = topk_accuracy(ctrl, "Top3")
top3_treat = topk_accuracy(treat, "Top3")
top5_ctrl = topk_accuracy(ctrl, "Top5")
top5_treat = topk_accuracy(treat, "Top5")

cats = ["Top-1", "Top-3", "Top-5"]
ctrl_vals = [top1_ctrl, top3_ctrl, top5_ctrl]
treat_vals = [top1_treat, top3_treat, top5_treat]

x = np.arange(len(cats))
w = 0.4

plt.figure(figsize=(6,4))
plt.bar(x - w/2, ctrl_vals, width=w, label="Control")
plt.bar(x + w/2, treat_vals, width=w, label="Treatment")
plt.xticks(x, cats)
plt.ylim(0, 1.0)
plt.ylabel("Accuracy / Hit rate")
plt.title(f"Top-k performance (source: {MODEL})")
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUTPUT}/{MODEL}_topk_summary_control_vs_treat.pdf", bbox_inches="tight")
plt.close()
