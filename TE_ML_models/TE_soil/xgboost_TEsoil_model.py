#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pickle
import argparse
import re

from xgboost import XGBClassifier

from sklearn.model_selection import GridSearchCV, LeaveOneGroupOut, cross_val_score
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    confusion_matrix,
    classification_report,
)

from skbio.stats.composition import clr

# Set up args
parser = argparse.ArgumentParser()
parser.add_argument("--model", choices=["control", "benzonase"], required=True,
    help="Which source treatment to train on (and which output folder/model name to use).")
parser.add_argument("--mode", choices=["train", "predict"], default="train",
                    help="train = fit and save model; predict = load saved model and run sinks/plots only")
parser.add_argument("--model_path", default=None,
                    help="Path to a saved .pkl model. If not set, uses default based on MODEL/OUTPUT.")
parser.add_argument("--taxonomy_tsv", default=None,
                    help="Optional QIIME2-style taxonomy.tsv with columns: 'Feature ID', 'Taxon', 'Confidence'. Adds Taxon to feature importance output.")
parser.add_argument("--run_permtest", action="store_true",
                    help="Run within-day permutation test on source LOGO CV.")
parser.add_argument("--threads", type=int, default=1,
                    help="Set number of threads to use.")
args = parser.parse_args()

# Set variables
MODEL = args.model
OUTPUT = f"{MODEL}_output"
default_model_path = f"{OUTPUT}/{MODEL}_best_estimator.pkl"
model_path = args.model_path if args.model_path else default_model_path
if args.threads < 1:
    raise ValueError("--threads must be >= 1")
THREADS = args.threads

os.makedirs(OUTPUT, exist_ok=True)

# Define path for the log file
log_path = f"{OUTPUT}/{MODEL}_log.txt"
logfile = open(log_path, "w")

# 0. Previous model bundle load if given
bundle = None
if args.mode == "predict":
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model bundle not found: {model_path}")
    with open(model_path, "rb") as f:
        bundle = pickle.load(f)

    best_model = bundle["model"]
    best_params = bundle.get("best_params", None)
    best_cv_score = bundle.get("best_cv_score", None)
    pseudocount_bundle = bundle.get("pseudocount", None)
    label_encoder_bundle = bundle.get("label_encoder", None)

    print(f"Loaded bundle from {model_path}")
    logfile.write(f"Loaded bundle from {model_path}\n")

#  1. Load feature table and metadata 
feature_table = pd.read_csv("TE_soil_table.tsv", sep="\t", index_col=0)

# If it looks like: rows = ASVs, columns = SampleIDs then transpose
if feature_table.columns.str.contains("NIJ-").any():
    feature_table = feature_table.transpose()

metadata = pd.read_csv("/home/zburcham/nij_analysis/ML_models/NIJ_metadata_model.txt", sep="\t")
metadata = metadata.set_index("SampleID")

valid_treatments = set(metadata["treatment"].astype(str).unique())
if MODEL not in valid_treatments:
    raise ValueError(f"--model {MODEL} not found in metadata['treatment']. Options: {sorted(valid_treatments)}")

shared_ids = metadata.index.intersection(feature_table.index)
shared_ids = shared_ids.sort_values()
metadata = metadata.loc[shared_ids]
feature_table = feature_table.loc[shared_ids]

assert metadata.shape[0] == feature_table.shape[0], "Mismatch in sample count"
assert all(metadata.index == feature_table.index), "Sample order mismatch"

print("Step 1: Loaded feature table and metadata.")

tax_lookup = None
if args.taxonomy_tsv:
    tax_df = pd.read_csv(args.taxonomy_tsv, sep="\t")
    # Normalize column names just in case
    cols = {c.strip(): c for c in tax_df.columns}
    if "Feature ID" not in cols or "Taxon" not in cols:
        raise ValueError(
            f"taxonomy_tsv must contain columns 'Feature ID' and 'Taxon'. Found: {list(tax_df.columns)}"
        )

    tax_df = tax_df.rename(columns={cols["Feature ID"]: "Feature ID", cols["Taxon"]: "Taxon"})
    tax_lookup = dict(zip(tax_df["Feature ID"].astype(str), tax_df["Taxon"].astype(str)))

    print(f"Loaded taxonomy table: {args.taxonomy_tsv} ({len(tax_lookup)} features)")
    logfile.write(f"Step 1.5:Loaded taxonomy table: {args.taxonomy_tsv} ({len(tax_lookup)} features)\n")

#  2. Split into source and sink before transformation 
source_mask = (metadata["mlm_1"] == "source") & (metadata["treatment"] == f"{MODEL}")
sink_mask_treatment = (metadata["mlm_1"] == "sink") & (metadata["treatment"] == "benzonase")
sink_mask_control = (metadata["mlm_1"] == "sink") & (metadata["treatment"] == "control")

X_source_raw = feature_table.loc[source_mask]
X_sink_raw_treatment = feature_table.loc[sink_mask_treatment]
X_sink_raw_control = feature_table.loc[sink_mask_control]

print("Step 2: Split into source and sink.")

assert not X_source_raw.empty, "Source (soil) feature table is empty!"
assert not X_sink_raw_treatment.empty, "Sink treatment feature table is empty!"
assert not X_sink_raw_control.empty, "Sink control feature table is empty!"

#  3. Pseudocount learned from TRAINING (source) only, then CLR
if args.mode == "predict":
    if pseudocount_bundle is None:
        raise ValueError("Bundle is missing 'pseudocount'. Re-train and save bundle.")
    pseudocount = float(pseudocount_bundle)
else:
    source_vals = X_source_raw.to_numpy()
    nonzero = source_vals[source_vals > 0]
    if nonzero.size == 0:
        raise ValueError("Source data has no non-zero entries; cannot compute pseudocount for CLR.")
    pseudocount = 0.5 * nonzero.min()

print(f"Step 3: Using pseudocount = {pseudocount}")
logfile.write(f"Pseudocount = {pseudocount}\n")

X_source_pc = X_source_raw.replace(0, pseudocount)
X_sink_pc_treatment = X_sink_raw_treatment.replace(0, pseudocount)
X_sink_pc_control = X_sink_raw_control.replace(0, pseudocount)

X_source_clr = clr(X_source_pc.values)
X_sink_clr_treatment = clr(X_sink_pc_treatment.values)
X_sink_clr_control = clr(X_sink_pc_control.values)

X_source = pd.DataFrame(X_source_clr, index=X_source_raw.index, columns=X_source_raw.columns)
X_sink_treatment = pd.DataFrame(X_sink_clr_treatment, index=X_sink_raw_treatment.index, columns=X_sink_raw_treatment.columns)
X_sink_control = pd.DataFrame(X_sink_clr_control, index=X_sink_raw_control.index, columns=X_sink_raw_control.columns)

print("Step 3: Applied CLR transformation.")

#  4. Labels and groups
y_source = metadata.loc[source_mask, "sample_location"].astype(str)
y_sink_treatment = metadata.loc[sink_mask_treatment, "sample_location"].astype(str)
y_sink_control = metadata.loc[sink_mask_control, "sample_location"].astype(str)

if args.mode == "predict":
    if label_encoder_bundle is None:
        raise ValueError("Bundle is missing 'label_encoder'. Re-train and save bundle.")
    label_encoder = label_encoder_bundle

    # sanity check: make sure current labels are known to the encoder
    unknown = set(pd.concat([y_source, y_sink_treatment, y_sink_control]).unique()) - set(label_encoder.classes_)
    if unknown:
        raise ValueError(f"Found labels not in bundle label_encoder: {sorted(unknown)}")

    y_source_encoded = label_encoder.transform(y_source)
    y_sink_encoded_treatment = label_encoder.transform(y_sink_treatment)
    y_sink_encoded_control = label_encoder.transform(y_sink_control)
else:
    label_encoder = LabelEncoder()
    y_source_encoded = label_encoder.fit_transform(y_source)
    y_sink_encoded_treatment = label_encoder.transform(y_sink_treatment)
    y_sink_encoded_control = label_encoder.transform(y_sink_control)

classes = label_encoder.classes_

print(f"Step 4: Labels encoded. Classes = {list(classes)}")

# Group for LeaveOneGroupOut (source only)
groups = metadata.loc[source_mask, "sampling_day"]
logo = LeaveOneGroupOut()

#  5. Model + hyperparams (tune for TOP-1 correctness) 
# For top-1 correctness, tune on balanced_accuracy (robust if classes are uneven).
xgb = XGBClassifier(
    objective="multi:softprob",
    eval_metric="mlogloss",
    random_state=42,
    n_jobs=1,
)

param_grid = {
    "n_estimators": [500],
    "learning_rate": [0.05, 0.1],
    "max_depth": [3, 6],
    "subsample": [0.8],
    "colsample_bytree": [0.8],
    "min_child_weight": [1, 3],
    "gamma": [0.0, 0.1],
}

grid_search = GridSearchCV(
    estimator=xgb,
    param_grid=param_grid,
    cv=logo.split(X_source, y_source_encoded, groups),
    scoring="balanced_accuracy",
    n_jobs=THREADS,
    refit=True,
)

print("Step 5: Configured GridSearchCV (scoring=balanced_accuracy).")

#  6. Fit on source samples only 
if args.mode == "train":
    print("Step 6: Fitting GridSearchCV...")
    grid_search.fit(X_source, y_source_encoded)

    best_model = grid_search.best_estimator_
    best_params = grid_search.best_params_
    best_cv_score = grid_search.best_score_

    print("Best Params:", best_params)
    print(f"Best CV balanced_accuracy (source LOGO): {best_cv_score:.3f}")

    logfile.write(f"Best Params: {best_params}\n")
    logfile.write(f"Best CV balanced_accuracy (source LOGO): {best_cv_score:.6f}\n")

    bundle = {
        "model": best_model,
        "best_params": best_params,
        "best_cv_score": float(best_cv_score),
        "label_encoder": label_encoder,
        "pseudocount": float(pseudocount),
        "classes": list(label_encoder.classes_),
    }
    with open(model_path, "wb") as f:
        pickle.dump(bundle, f)
    print(f"Model saved to {model_path}")
    logfile.write(f"Model saved to {model_path}\n")

else:
    # Predict mode sanity checks (bundle should have been loaded earlier)
    if bundle is None:
        raise RuntimeError(
            "Predict mode requires a loaded bundle, but bundle is None. "
            "Did you load it before preprocessing? Check your early bundle-load block."
        )
    if "model" not in bundle or best_model is None:
        raise RuntimeError(
            "Predict mode bundle is missing the trained model (key 'model') or best_model is None."
        )
    if bundle.get("pseudocount", None) is None:
        raise RuntimeError(
            "Predict mode bundle is missing 'pseudocount'. Re-train and save a bundle including pseudocount."
        )
    if bundle.get("label_encoder", None) is None:
        raise RuntimeError(
            "Predict mode bundle is missing 'label_encoder'. Re-train and save a bundle including label_encoder."
        )

    print("Step 6: Predict mode (bundle loaded; model/preproc confirmed).")
    logfile.write("Step 6: Predict mode (bundle loaded; model/preproc confirmed).\n")

    if best_params is not None:
        print("Loaded best_params:", best_params)
        logfile.write(f"Loaded best_params: {best_params}\n")

# Feature importances
importances = best_model.feature_importances_
feature_names = X_source.columns.astype(str)

importance_df = pd.DataFrame(
    {"Feature": feature_names, "Importance": importances}
).sort_values(by="Importance", ascending=False)

# Add taxonomy if provided
if tax_lookup is not None:
    importance_df["Taxon"] = importance_df["Feature"].map(tax_lookup)
    missing = importance_df["Taxon"].isna().sum()
    if missing > 0:
        print(f"Warning: {missing} features missing taxonomy annotation (will be NaN).")
        logfile.write(f"Warning: {missing} features missing taxonomy annotation (will be NaN).\n")

importance_df.to_csv(f"{OUTPUT}/{MODEL}_feature_importance.csv", index=False)

print("Model fitting/loading complete.")

#  7. Predict on sink samples 
print("Step 7: Predicting on sink samples...")

y_pred_treatment = best_model.predict(X_sink_treatment)
decoded_preds_treatment = label_encoder.inverse_transform(y_pred_treatment)

y_pred_control = best_model.predict(X_sink_control)
decoded_preds_control = label_encoder.inverse_transform(y_pred_control)

predictions_df_treatment = pd.DataFrame(
    {"SampleID": X_sink_treatment.index, "True_Label": y_sink_treatment.values, "Predicted_Label": decoded_preds_treatment}
)
predictions_df_treatment.to_csv(f"{OUTPUT}/source_{MODEL}_sink_treatment_predictions.csv", index=False)

predictions_df_control = pd.DataFrame(
    {"SampleID": X_sink_control.index, "True_Label": y_sink_control.values, "Predicted_Label": decoded_preds_control}
)
predictions_df_control.to_csv(f"{OUTPUT}/source_{MODEL}_sink_control_predictions.csv", index=False)

print("Predicted on sink samples.")

#  8. Evaluate TOP-1 correctness 
print("Step 8: Evaluating TOP-1 location calling...")

def eval_and_log(name, y_true_int, y_pred_int, y_true_str, y_pred_str):
    acc = accuracy_score(y_true_int, y_pred_int)
    bacc = balanced_accuracy_score(y_true_int, y_pred_int)
    f1m = f1_score(y_true_int, y_pred_int, average="macro")

    report = classification_report(
        y_true_str,
        y_pred_str,
        labels=classes,
        target_names=classes,
        digits=3,
        zero_division=0,
    )

    print(f"\n{name}")
    print(f"  Accuracy:          {acc:.3f}")
    print(f"  Balanced accuracy: {bacc:.3f}")
    print(f"  Macro-F1:          {f1m:.3f}")
    print(report)

    logfile.write(f"\n{name}\n")
    logfile.write(f"Accuracy: {acc:.6f}\n")
    logfile.write(f"Balanced accuracy: {bacc:.6f}\n")
    logfile.write(f"Macro-F1: {f1m:.6f}\n")
    logfile.write(report + "\n")

    return acc, bacc, f1m

# Treatment sink
_ = eval_and_log(
    "Treatment sink",
    y_sink_encoded_treatment,
    y_pred_treatment,
    y_sink_treatment.values,
    decoded_preds_treatment,
)

# Control sink
_ = eval_and_log(
    "Control sink",
    y_sink_encoded_control,
    y_pred_control,
    y_sink_control.values,
    decoded_preds_control,
)

print("Evaluated predictions.")

# --- Confusion matrices (normalized by true label) ---
print("Plotting confusion matrices...")

cm_treatment = confusion_matrix(
    y_sink_treatment.values, decoded_preds_treatment,
    labels=classes, normalize="true"
)
cm_control = confusion_matrix(
    y_sink_control.values, decoded_preds_control,
    labels=classes, normalize="true"
)

display_labels = [c.replace("_", " ").title() for c in classes]

plt.figure(figsize=(8, 6))
sns.heatmap(cm_treatment, annot=True, fmt=".2f", cmap="Blues", xticklabels=display_labels, yticklabels=display_labels)
plt.ylabel("True label")
plt.xlabel("Predicted label")
plt.title("Normalized Confusion Matrix (source: , sink: treatment)")
plt.tight_layout()
plt.savefig(f"{OUTPUT}/{MODEL}_on_treatment_confusion_matrix_normalized.pdf", bbox_inches="tight")
plt.close()

plt.figure(figsize=(8, 6))
sns.heatmap(cm_control, annot=True, fmt=".2f", cmap="Blues", xticklabels=display_labels, yticklabels=display_labels)
plt.ylabel("True label")
plt.xlabel("Predicted label")
plt.title("Normalized Confusion Matrix (source: , sink: control)")
plt.tight_layout()
plt.savefig(f"{OUTPUT}/{MODEL}_on_control_confusion_matrix_normalized.pdf", bbox_inches="tight")
plt.close()

#  10. Plot important features 
top_n = 10
top_features = importance_df.head(top_n)

plt.figure(figsize=(8, 6))
plt.barh(top_features["Feature"], top_features["Importance"])
plt.gca().invert_yaxis()
plt.xlabel("Feature Importance")
plt.title(f"Top {top_n} Features — {MODEL}")
plt.tight_layout()
plt.savefig(f"{OUTPUT}/{MODEL}_top{top_n}_features.pdf", bbox_inches="tight")
plt.close()

print("Generated figures.")

####################################################################################################################################

# Permutation test (null model test)
if args.run_permtest:
    def permute_within_groups(y, groups, rng):
        """Permute labels within each group (sampling_day)."""
        y_perm = y.copy()
        groups = np.asarray(groups)
        for g in np.unique(groups):
            idx = np.where(groups == g)[0]
            y_perm[idx] = rng.permutation(y_perm[idx])
        return y_perm

    def run_logo_cv_balacc(X, y, groups, xgb_params, n_jobs=1):
        """Compute LOGO CV balanced accuracy for fixed params."""
        logo = LeaveOneGroupOut()
        clf = XGBClassifier(
            objective="multi:softprob",
            eval_metric="mlogloss",
            random_state=42,
            n_jobs=1, 
            **xgb_params
        )
        scores = cross_val_score(
            clf,
            X,
            y,
            cv=logo.split(X, y, groups),
            scoring="balanced_accuracy",
            n_jobs=n_jobs
        )
        return scores

    # 1) Observed CV score (real labels) using fixed best_params
    obs_scores = run_logo_cv_balacc(X_source, y_source_encoded, groups, best_params, n_jobs=THREADS)
    obs_mean = obs_scores.mean()

    print(f"\nPermutation test setup:")
    print(f"Observed LOGO CV balanced_accuracy mean = {obs_mean:.6f}")
    logfile.write(f"\nObserved LOGO CV balanced_accuracy mean = {obs_mean:.6f}\n")

    # 2) Permutations
    N_PERM = 999
    rng = np.random.default_rng(2026)

    perm_means = np.empty(N_PERM, dtype=float)

    for i in range(N_PERM):
        y_perm = permute_within_groups(y_source_encoded, groups, rng)
        perm_scores = run_logo_cv_balacc(X_source, y_perm, groups, best_params, n_jobs=THREADS)
        perm_means[i] = perm_scores.mean()
        if (i + 1) % 50 == 0:
            print(f"  completed {i+1}/{N_PERM} permutations")

    # 3) One-sided p-value: P(null >= observed)
    # Add +1 smoothing so p is never exactly 0
    p_value = (1 + np.sum(perm_means >= obs_mean)) / (N_PERM + 1)

    print(f"\nPermutation test results (within-day label shuffle):")
    print(f"Observed mean bal_acc = {obs_mean:.6f}")
    print(f"Null mean bal_acc     = {perm_means.mean():.6f}")
    print(f"Null max bal_acc      = {perm_means.max():.6f}")
    print(f"p-value (>= obs)      = {p_value:.6g}")

    logfile.write("\nPermutation test results (within-day label shuffle):\n")
    logfile.write(f"Observed mean bal_acc = {obs_mean:.6f}\n")
    logfile.write(f"Null mean bal_acc     = {perm_means.mean():.6f}\n")
    logfile.write(f"Null max bal_acc      = {perm_means.max():.6f}\n")
    logfile.write(f"p-value (>= obs)      = {p_value:.6g}\n")

    # 4) Save results
    perm_df = pd.DataFrame({"perm_mean_balacc": perm_means})
    perm_out = f"{OUTPUT}/{MODEL}_permtest_within_day_balacc.csv"
    perm_df.to_csv(perm_out, index=False)
    print(f"Saved permutation null distribution to {perm_out}")
    logfile.write(f"Saved permutation null distribution to {perm_out}\n")

    # 5) Plot histogram
    plt.figure(figsize=(7, 5))
    plt.hist(perm_means, bins=30)
    plt.axvline(obs_mean, linewidth=2)
    plt.xlabel("Mean LOGO CV balanced accuracy (permuted labels)")
    plt.ylabel("Count")
    plt.title("Permutation test null distribution (within sampling_day)")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT}/{MODEL}_permtest_within_day_balacc_hist.pdf", bbox_inches="tight")
    plt.close()

logfile.close()
