#!/usr/bin/env python3

import argparse
import os
import pandas as pd
from math import comb

def topk_hit(true_label: str, topk_str: str, sep: str = ";") -> int:
    if pd.isna(topk_str):
        return 0
    items = [x.strip() for x in str(topk_str).split(sep) if x.strip() != ""]
    return int(str(true_label).strip() in items)

def mcnemar_exact_pvalue(b: int, c: int) -> float:
    """Exact two-sided McNemar p-value via binomial test on discordant pairs."""
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    cdf = sum(comb(n, i) for i in range(0, k + 1)) / (2 ** n)
    return min(1.0, 2 * cdf)


def paired_test_from_predictions(
    control_csv: str,
    treatment_csv: str,
    out_dir: str,
    paired_csv_path: str,
    log_path: str,
    strip_regex: str,
):
    os.makedirs(out_dir, exist_ok=True)

    df_c = pd.read_csv(control_csv)
    df_t = pd.read_csv(treatment_csv)

    # Sanity: strip whitespace
    df_c["SampleID"] = df_c["SampleID"].astype(str).str.strip()
    df_t["SampleID"] = df_t["SampleID"].astype(str).str.strip()

    # Pair id: remove trailing A/B (customizable)
    df_c["pair_id"] = df_c["SampleID"].str.replace(strip_regex, "", regex=True)
    df_t["pair_id"] = df_t["SampleID"].str.replace(strip_regex, "", regex=True)

    # Correctness
    if "Top3" not in df_c.columns or "Top3" not in df_t.columns:
        raise ValueError("Top3 column not found in one or both CSVs. Re-run prediction script to write Top3.")

    df_c["correct1_control"] = (df_c["True_Label"] == df_c["Predicted_Label"]).astype(int)
    df_t["correct1_treat"] = (df_t["True_Label"] == df_t["Predicted_Label"]).astype(int)

    df_c["correct3_control"] = [topk_hit(t, k) for t, k in zip(df_c["True_Label"], df_c["Top3"])]
    df_t["correct3_treat"] = [topk_hit(t, k) for t, k in zip(df_t["True_Label"], df_t["Top3"])]

    # using top-3 for test
    df_c["correct_control"] = df_c["correct3_control"]
    df_t["correct_treat"] = df_t["correct3_treat"]

    # Merge pairs
    pairs = df_c[["pair_id", "SampleID", "True_Label", "Predicted_Label", "correct_control"]].merge(
        df_t[["pair_id", "SampleID", "True_Label", "Predicted_Label", "correct_treat"]],
        on="pair_id",
        how="inner",
        suffixes=("_control", "_treat"),
    )

    n_pairs = len(pairs)
    if n_pairs == 0:
        raise ValueError(
            "No pairs matched. Check SampleID naming and the --strip_regex rule. "
            f"Current strip_regex={strip_regex!r}"
        )

    # Check for mismatched True_Label
    mism_true = int((pairs["True_Label_control"] != pairs["True_Label_treat"]).sum())

    # Discordant counts
    b = int(((pairs["correct_control"] == 1) & (pairs["correct_treat"] == 0)).sum())  # control-only correct
    c = int(((pairs["correct_control"] == 0) & (pairs["correct_treat"] == 1)).sum())  # treatment-only correct

    acc_control = float(pairs["correct_control"].mean())
    acc_treat = float(pairs["correct_treat"].mean())
    delta = acc_treat - acc_control

    p_exact = mcnemar_exact_pvalue(b, c)

    # Write paired table
    pairs.to_csv(paired_csv_path, index=False)

    # Log + print
    def logprint(logf, *args):
        line = " ".join(str(a) for a in args)
        print(line)
        logf.write(line + "\n")
        logf.flush()

    with open(log_path, "a") as logf:
        logprint(logf, "\n=== Paired McNemar Test ===")
        logprint(logf, "Control file:", control_csv)
        logprint(logf, "Treatment file:", treatment_csv)
        logprint(logf, "Paired n =", n_pairs)
        if mism_true:
            logprint(logf, f"WARNING: {mism_true} pairs have different True_Label between control and treatment.")
        logprint(logf, f"Accuracy control (B):   {acc_control:.3f}")
        logprint(logf, f"Accuracy treatment (A): {acc_treat:.3f}")
        logprint(logf, f"Δ accuracy (A - B):     {delta:.3f}")
        logprint(logf, "Discordant pairs:")
        logprint(logf, "  b = control correct, treatment wrong =", b)
        logprint(logf, "  c = control wrong, treatment correct =", c)
        logprint(logf, f"McNemar exact p-value (two-sided): {p_exact:.6g}")
        logprint(logf, "Wrote paired table:", paired_csv_path)

    return pairs, {
        "n_pairs": n_pairs,
        "acc_control": acc_control,
        "acc_treat": acc_treat,
        "delta": delta,
        "b": b,
        "c": c,
        "p": p_exact,
        "mism_true": mism_true,
        "paired_csv": paired_csv_path,
        "log": log_path,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Paired accuracy comparison (control vs treatment) using an exact McNemar test."
    )
    parser.add_argument("--control_csv", required=True, help="CSV with SampleID, True_Label, Predicted_Label (control).")
    parser.add_argument("--treatment_csv", required=True, help="CSV with SampleID, True_Label, Predicted_Label (treatment).")
    parser.add_argument(
        "--out_dir",
        default=None,
        help="Output directory for log and paired CSV. Default: directory of --control_csv.",
    )
    parser.add_argument(
        "--log",
        default=None,
        help="Log file path. Default: <out_dir>/paired_test.log",
    )
    parser.add_argument(
        "--paired_csv",
        default=None,
        help="Output paired table path. Default: <out_dir>/paired_control_vs_treatment_predictions.csv",
    )
    parser.add_argument(
        "--strip_regex",
        default=r"[AB]$",
        help=r"Regex to strip from SampleID to form pair_id. Default: '[AB]$' (removes trailing A/B).",
    )

    args = parser.parse_args()

    out_dir = args.out_dir if args.out_dir else os.path.dirname(args.control_csv) or "."
    log_path = args.log if args.log else os.path.join(out_dir, "paired_test.log")
    paired_csv_path = args.paired_csv if args.paired_csv else os.path.join(
        out_dir, "paired_control_vs_treatment_predictions.csv"
    )

    paired_test_from_predictions(
        control_csv=args.control_csv,
        treatment_csv=args.treatment_csv,
        out_dir=out_dir,
        paired_csv_path=paired_csv_path,
        log_path=log_path,
        strip_regex=args.strip_regex,
    )


if __name__ == "__main__":
    main()
