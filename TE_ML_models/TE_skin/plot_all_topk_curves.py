#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt

# Point these at your four CSVs
files = {
    "Source: Control, Sink: Control":   "control_output/control_topk_curve_control.csv",
    "Source: Control, Sink: Treatment": "control_output/control_topk_curve_treatment.csv",
    "Source: Treatment, Sink: Control":       "benzonase_output/benzonase_topk_curve_control.csv",
    "Source: Treatment, Sink: Treatment":     "benzonase_output/benzonase_topk_curve_treatment.csv",
}

style = {
    "Source: Control, Sink: Control":     dict(color="tab:blue",   linestyle="-"),
    "Source: Control, Sink: Treatment":   dict(color="tab:blue",   linestyle="--"),
    "Source: Treatment, Sink: Control":   dict(color="tab:orange", linestyle="-"),
    "Source: Treatment, Sink: Treatment": dict(color="tab:orange", linestyle="--"),
}

plt.figure(figsize=(7, 4))

for label, path in files.items():
    df = pd.read_csv(path)
    plt.plot(
        df["k"], df["hit_rate"],
        marker="o", markersize=3,
        label=label,
        **style[label]
    )

plt.ylim(0, 1.0)
plt.xlabel("k (Top-k)")
plt.ylabel("Top-k hit rate")
plt.title("Top-k curves across models and conditions")
plt.xticks(rotation=90)
plt.legend(fontsize=8, labelspacing=0.3, handlelength=2.5, handletextpad=0.6, borderpad=0.3)
plt.tight_layout()
plt.savefig("all_models_topk_curves.pdf", bbox_inches="tight")
plt.close()
