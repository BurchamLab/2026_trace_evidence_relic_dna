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

from sklearn.model_selection import GroupKFold, GridSearchCV, StratifiedKFold, StratifiedGroupKFold, LeaveOneGroupOut, cross_val_score
from sklearn.preprocessing import LabelBinarizer, label_binarize, LabelEncoder
from sklearn.metrics import accuracy_score, roc_auc_score, ConfusionMatrixDisplay, confusion_matrix, RocCurveDisplay, roc_curve, auc, balanced_accuracy_score, f1_score, classification_report, recall_score
from skbio.stats.composition import clr, multi_replace
from collections import defaultdict


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
                    help="Run global permutation test on source SKG CV.")
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
    label_encoder = bundle.get("label_encoder", None)
    if label_encoder is None:
        raise ValueError("Bundle missing 'label_encoder'. Re-train to create a compatible bundle.")

    print(f"Loaded bundle from {model_path}")
    logfile.write(f"Loaded bundle from {model_path}\n")

#  1. Load feature table and metadata 
feature_table = pd.read_csv('TE_touch_table.tsv', sep='\t', index_col=0)

# If it looks like: rows = ASVs, columns = SampleIDs then transpose
if feature_table.columns.str.contains("NIJ-").any():
    feature_table = feature_table.transpose()

metadata = pd.read_csv("/home/zburcham/nij_analysis/ML_models/NIJ_metadata_model.txt", sep="\t")
metadata = metadata.set_index('SampleID')
metadata["pair_id"] = metadata.index.to_series().str.replace(r"[AB]$", "", regex=True)

valid_treatments = set(metadata["treatment"].astype(str).str.strip().str.lower().unique())
if MODEL.lower() not in valid_treatments:
    raise ValueError(f"--model {MODEL} not found in metadata['treatment']. Options: {sorted(valid_treatments)}")

metadata.index = metadata.index.astype(str).str.strip()
feature_table.index = feature_table.index.astype(str).str.strip()
feature_table.columns = feature_table.columns.astype(str).str.strip()

shared_ids = metadata.index.intersection(feature_table.index)
shared_ids = shared_ids.sort_values()
metadata = metadata.loc[shared_ids]
feature_table = feature_table.loc[shared_ids]

assert metadata.shape[0] == feature_table.shape[0], "Mismatch in sample count"
assert all(metadata.index == feature_table.index), "Sample order mismatch"

if len(shared_ids) == 0:
    raise ValueError("No shared SampleIDs between metadata and feature table after alignment.")

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

# Normalize text columns (prevents silent mask failures)
metadata["treatment"] = metadata["treatment"].astype(str).str.strip().str.lower()
metadata["sample_type"] = metadata["sample_type"].astype(str).str.strip().str.lower()
metadata["palm"] = metadata["palm"].astype(str).str.strip().str.lower()

# Create the class label = person + palm and just person
metadata["person_palm"] = metadata["anon_donor_ID"].astype(str) + "_" + metadata["palm"]
metadata["person_label"] = metadata["anon_donor_ID"].astype(str)

#  2. Split into source and sink before transformation 
source_mask = (metadata["sample_type"] == "skin") & (metadata["treatment"] == f"{MODEL}")
sink_mask_control = (metadata["sample_type"] == "pvc") & (metadata["treatment"] == "control")
sink_mask_treatment = (metadata["sample_type"] == "pvc") & (metadata["treatment"] == "benzonase")

X_source_raw = feature_table.loc[source_mask]
X_sink_raw_control = feature_table.loc[sink_mask_control]
X_sink_raw_treatment = feature_table.loc[sink_mask_treatment]

y_source = metadata.loc[source_mask, "person_label"].astype(str)
y_sink_control = metadata.loc[sink_mask_control, "person_label"].astype(str)
y_sink_treatment = metadata.loc[sink_mask_treatment, "person_label"].astype(str)

assert not X_source_raw.empty, "No source samples found after filtering."
assert not X_sink_raw_control.empty, "No sink PVC control samples found."
assert not X_sink_raw_treatment.empty, "No sink PVC benzonase samples found."

print("Counts:",
      "source(skin)=", X_source_raw.shape[0],
      "sink(pvc control)=", X_sink_raw_control.shape[0],
      "sink(pvc benzonase)=", X_sink_raw_treatment.shape[0])
      
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

# 4. Labels (encode for XGBoost; keep human-readable class names for reporting)
if args.mode == "train":
    label_encoder = LabelEncoder()
    y_source_enc = label_encoder.fit_transform(y_source.values)
else:
    # predict mode: encoder loaded from bundle above
    y_source_enc = label_encoder.transform(y_source.values)

classes = label_encoder.classes_  # donor IDs in encoder order
print(f"Step 4: Encoded labels. n_classes = {len(classes)}")
logfile.write(f"Step 4: Encoded labels. n_classes = {len(classes)}\n")

# 5. Set up model + grouped CV splits + GridSearchCV
def donor_stratified_group_splits(y, groups, n_splits=2, seed=42):
    """
    Build CV splits at the GROUP (pair_id) level while stratifying by class label (y).

    Goal:
      - Keep each group (pair_id) intact to prevent leakage.
      - Distribute each class's groups across folds so that (when feasible) every class
        is present in every training fold ("closed-set" grouped CV).

    Feasibility:
      - Requires each class to have >= n_splits unique groups. If not, raises ValueError.

    Assumption:
      - Each group (pair_id) contains samples from exactly one class.
    """
    rng = np.random.default_rng(seed)
    y = np.asarray(y)
    groups = np.asarray(groups)

    # group -> donor (must be unique)
    group_to_label = {}
    for g in np.unique(groups):
        idx = np.where(groups == g)[0]
        donors = np.unique(y[idx])
        if len(donors) != 1:
            raise ValueError(f"pair_id '{g}' contains multiple donors: {donors}")
        group_to_label[g] = donors[0]

    label_to_groups = defaultdict(list)
    for g, d in group_to_label.items():
        label_to_groups[d].append(g)

    min_g = min(len(gs) for gs in label_to_groups.values())
    if n_splits > min_g:
        raise ValueError(
            f"n_splits={n_splits} too large for grouped CV: smallest class has only {min_g} pair_id groups."
        )

    # assign groups to folds donor-by-donor (round-robin)
    fold_groups = [set() for _ in range(n_splits)]
    for d, gs in label_to_groups.items():
        gs = np.array(gs)
        rng.shuffle(gs)
        for i, g in enumerate(gs):
            fold_groups[i % n_splits].add(g)

    splits = []
    for k in range(n_splits):
        test_g = fold_groups[k]
        test_idx = np.where(np.isin(groups, list(test_g)))[0]
        train_idx = np.where(~np.isin(groups, list(test_g)))[0]
        splits.append((train_idx, test_idx))
    return splits, min_g

groups = metadata.loc[source_mask, "pair_id"].values
y = y_source_enc

# Feasible n_splits based on min unique groups per donor
tmp = pd.DataFrame({"y": y, "g": groups})
min_groups_per_class = int(tmp.groupby("y")["g"].nunique().min())
# Strongly recommend starting with 2
n_splits = min(2, min_groups_per_class)

if n_splits < 2:
    raise ValueError(
        f"Not enough pair_id groups per donor for grouped CV. min_groups_per_class={min_groups_per_class}."
    )

cv_splits, min_g = donor_stratified_group_splits(y, groups, n_splits=n_splits, seed=42)
assert min_g == min_groups_per_class

all_labels = set(y)
for i, (tr, te) in enumerate(cv_splits):
    missing = all_labels - set(y[tr])
    if missing:
        raise ValueError(f"Fold {i} missing classes in TRAIN: {sorted(missing)}")

print(f"CV: using donor-stratified grouped splits with n_splits={n_splits} (min groups per donor = {min_groups_per_class})")
logfile.write(f"CV: using donor-stratified grouped splits with n_splits={n_splits} (min groups per donor = {min_groups_per_class})\n")

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
    cv=cv_splits,
    scoring="balanced_accuracy",
    n_jobs=THREADS,
    refit=True,
)

print("Step 5: Configured XGBoost.")

#  6. Fit on source samples only 
if args.mode == "train":
    print("Step 6: Fitting GridSearchCV...")
    grid_search.fit(X_source, y_source_enc)

    best_model = grid_search.best_estimator_
    best_params = grid_search.best_params_
    best_cv_score = grid_search.best_score_

    print("Best Params:", best_params)
    print(f"Best CV balanced_accuracy (source SKG): {best_cv_score:.3f}")

    logfile.write(f"Best Params: {best_params}\n")
    logfile.write(f"Best CV balanced_accuracy (source SKG): {best_cv_score:.6f}\n")

    bundle = {
        "model": best_model,
        "best_params": best_params,
        "best_cv_score": float(best_cv_score),
        "pseudocount": float(pseudocount),
        "label_encoder": label_encoder
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

# Top-1 predictions (encoded)
y_pred_treatment_enc = best_model.predict(X_sink_treatment)
y_pred_control_enc   = best_model.predict(X_sink_control)

# Probabilities for Top-k
y_proba_treatment = best_model.predict_proba(X_sink_treatment)
y_proba_control   = best_model.predict_proba(X_sink_control)

# Convert Top-1 to original string labels
y_pred_treatment = label_encoder.inverse_transform(y_pred_treatment_enc.astype(int))
y_pred_control   = label_encoder.inverse_transform(y_pred_control_enc.astype(int))

# ---- Top-k label helper ----
def top_k_labels(y_proba, model_classes, label_encoder, k=3):
    """
    y_proba: (n_samples, n_classes) probabilities aligned to model_classes
    model_classes: best_model.classes_ (encoded class IDs, in proba column order)
    Returns: array of shape (n_samples, k) of decoded string labels
    """
    model_classes = np.asarray(model_classes)
    topk_idx = np.argsort(y_proba, axis=1)[:, -k:][:, ::-1]  # descending
    topk_enc = model_classes[topk_idx]  # encoded class IDs
    # decode each row
    return np.array([label_encoder.inverse_transform(row.astype(int)) for row in topk_enc])


# proba columns align to these classes:
model_classes_enc = best_model.classes_

top2_labels_treat   = top_k_labels(y_proba_treatment, model_classes_enc, label_encoder, k=2)
top2_labels_control = top_k_labels(y_proba_control,   model_classes_enc, label_encoder, k=2)

top3_labels_treat   = top_k_labels(y_proba_treatment, model_classes_enc, label_encoder, k=3)
top3_labels_control = top_k_labels(y_proba_control,   model_classes_enc, label_encoder, k=3)

top5_labels_treat   = top_k_labels(y_proba_treatment, model_classes_enc, label_encoder, k=5)
top5_labels_control = top_k_labels(y_proba_control,   model_classes_enc, label_encoder, k=5)

# Build outputs
predictions_df_treatment = pd.DataFrame(
    {
        "SampleID": X_sink_treatment.index,
        "True_Label": y_sink_treatment.values,
        "Predicted_Label": y_pred_treatment,
        "Top2": [";".join(row) for row in top2_labels_treat],
        "Top3": [";".join(row) for row in top3_labels_treat],
        "Top5": [";".join(row) for row in top5_labels_treat],
    }
)

predictions_df_treatment.to_csv(f"{OUTPUT}/source_{MODEL}_sink_treatment_predictions.csv", index=False)

predictions_df_control = pd.DataFrame(
    {
        "SampleID": X_sink_control.index,
        "True_Label": y_sink_control.values,
        "Predicted_Label": y_pred_control,
        "Top2": [";".join(row) for row in top2_labels_control],
        "Top3": [";".join(row) for row in top3_labels_control],
        "Top5": [";".join(row) for row in top5_labels_control],
    }
)
predictions_df_control.to_csv(f"{OUTPUT}/source_{MODEL}_sink_control_predictions.csv", index=False)

print("Predicted on sink samples.")

#  8. Evaluate 
print("Step 8: Evaluating predictions....")

accuracy_treatment = accuracy_score(y_sink_treatment, y_pred_treatment)
accuracy_control = accuracy_score(y_sink_control, y_pred_control)

def macro_auc_ovr_strings(y_true, y_proba, classes):
    # y_proba columns are aligned to encoded class indices (0..K-1), which map to label_encoder.classes_
    Y = label_binarize(y_true, classes=classes)
    aucs = []
    for i in range(len(classes)):
        pos = Y[:, i].sum()
        if pos == 0 or pos == Y.shape[0]:
            continue
        aucs.append(roc_auc_score(Y[:, i], y_proba[:, i]))
    return float(np.mean(aucs)) if aucs else np.nan


# proba columns correspond to encoded class indices 0..K-1, which map to label_encoder.classes_
proba_classes = label_encoder.inverse_transform(model_classes_enc.astype(int))

macro_auc_treatment = macro_auc_ovr_strings(y_sink_treatment, y_proba_treatment, proba_classes)
macro_auc_control   = macro_auc_ovr_strings(y_sink_control,   y_proba_control,   proba_classes)

balacc_treatment = balanced_accuracy_score(y_sink_treatment, y_pred_treatment)
balacc_control   = balanced_accuracy_score(y_sink_control, y_pred_control)

macro_f1_treatment = f1_score(y_sink_treatment, y_pred_treatment, average="macro")
macro_f1_control   = f1_score(y_sink_control, y_pred_control, average="macro")

print(f"Accuracy on PVC benzonase: {accuracy_treatment:.3f}")
print(f"PVC benzonase Macro AUC (OvR): {macro_auc_treatment:.3f}")
print(f"PVC benzonase Balanced accuracy: {balacc_treatment:.3f} | Macro-F1: {macro_f1_treatment:.3f}")
print(f"Accuracy on PVC control: {accuracy_control:.3f}")
print(f"PVC control Macro AUC (OvR): {macro_auc_control:.3f}")
print(f"PVC control   Balanced accuracy: {balacc_control:.3f} | Macro-F1: {macro_f1_control:.3f}")

# check top 3 and top 5 accuracy
def top_k_accuracy_strings(y_true, y_proba, classes, k=3):
    """
    y_true: array-like of true class labels (strings)
    y_proba: (n_samples, n_classes) probabilities aligned to `classes`
    classes: array-like of class labels (strings) in proba column order
    """
    y_true = np.asarray(y_true, dtype=str)
    classes = np.asarray(classes, dtype=str)

    # indices of top-k probs per sample (largest k)
    topk_idx = np.argsort(y_proba, axis=1)[:, -k:]
    # map indices -> labels
    topk_labels = classes[topk_idx]
    # membership test per row
    hits = [y_true[i] in topk_labels[i] for i in range(len(y_true))]
    return float(np.mean(hits))

top2_treatment = top_k_accuracy_strings(y_sink_treatment, y_proba_treatment, proba_classes, k=2)
top2_control   = top_k_accuracy_strings(y_sink_control,   y_proba_control,   proba_classes, k=2)

print(f"PVC benzonase Top-2 accuracy: {top2_treatment:.3f}")
print(f"PVC control   Top-2 accuracy: {top2_control:.3f}")

logfile.write(f"PVC benzonase Top-3 accuracy: {top2_treatment:.3f}\n")
logfile.write(f"PVC control   Top-3 accuracy: {top2_control:.3f}\n")

top3_treatment = top_k_accuracy_strings(y_sink_treatment, y_proba_treatment, proba_classes, k=3)
top3_control   = top_k_accuracy_strings(y_sink_control,   y_proba_control,   proba_classes, k=3)

print(f"PVC benzonase Top-3 accuracy: {top3_treatment:.3f}")
print(f"PVC control   Top-3 accuracy: {top3_control:.3f}")

logfile.write(f"PVC benzonase Top-3 accuracy: {top3_treatment:.3f}\n")
logfile.write(f"PVC control   Top-3 accuracy: {top3_control:.3f}\n")

top5_treatment = top_k_accuracy_strings(y_sink_treatment, y_proba_treatment, proba_classes, k=5)
top5_control   = top_k_accuracy_strings(y_sink_control,   y_proba_control,   proba_classes, k=5)

print(f"PVC benzonase Top-5 accuracy: {top5_treatment:.3f}")
print(f"PVC control   Top-5 accuracy: {top5_control:.3f}")

logfile.write(f"PVC benzonase Top-5 accuracy: {top5_treatment:.3f}\n")
logfile.write(f"PVC control   Top-5 accuracy: {top5_control:.3f}\n")

# Write to log
logfile.write(f"Accuracy on PVC benzonase: {accuracy_treatment:.3f}\n")
logfile.write(f"PVC benzonase Macro AUC (OvR): {macro_auc_treatment:.3f}\n")
logfile.write(f"PVC benzonase Balanced accuracy: {balacc_treatment:.3f} | Macro-F1: {macro_f1_treatment:.3f}\n")
logfile.write(f"Accuracy on PVC control: {accuracy_control:.3f}\n")
logfile.write(f"PVC control Macro AUC (OvR): {macro_auc_control:.3f}\n")
logfile.write(f"PVC control   Balanced accuracy: {balacc_control:.3f} | Macro-F1: {macro_f1_control:.3f}\n")

# Adding reject option
def top1_top2_and_margin(y_proba):
    """
    Returns p1 (top1 prob), p2 (top2 prob), margin (p1-p2), and top1 index.
    """
    # argsort ascending; take last two
    idx_sorted = np.argsort(y_proba, axis=1)
    top1_idx = idx_sorted[:, -1]
    top2_idx = idx_sorted[:, -2]

    p1 = y_proba[np.arange(y_proba.shape[0]), top1_idx]
    p2 = y_proba[np.arange(y_proba.shape[0]), top2_idx]
    margin = p1 - p2
    return p1, p2, margin, top1_idx

def reject_metrics(y_true, y_pred, score, threshold):
    """
    Reject when score < threshold.
    Returns coverage, accuracy_kept, n_kept, n_total.
    """
    y_true = np.asarray(y_true, dtype=str)
    y_pred = np.asarray(y_pred, dtype=str)
    score = np.asarray(score, dtype=float)

    keep = score >= threshold
    n_total = len(y_true)
    n_kept = int(np.sum(keep))
    coverage = n_kept / n_total if n_total else np.nan

    acc_kept = np.mean(y_true[keep] == y_pred[keep]) if n_kept > 0 else np.nan
    return coverage, float(acc_kept), n_kept, n_total, keep

model_classes_enc = best_model.classes_
proba_classes = label_encoder.inverse_transform(model_classes_enc.astype(int))

# Treatment scores
p1_t, p2_t, margin_t, top1_idx_t = top1_top2_and_margin(y_proba_treatment)
top1_label_t = proba_classes[top1_idx_t]

# Control scores
p1_c, p2_c, margin_c, top1_idx_c = top1_top2_and_margin(y_proba_control)
top1_label_c = proba_classes[top1_idx_c]

print("\nReject option using margin = p_top1 - p_top2")
for thr in [0.00, 0.05, 0.10, 0.15, 0.20]:
    cov_t, acc_t, nk_t, nt_t, _ = reject_metrics(y_sink_treatment, top1_label_t, margin_t, thr)
    cov_c, acc_c, nk_c, nt_c, _ = reject_metrics(y_sink_control,   top1_label_c, margin_c, thr)

    print(f"  thr={thr:0.2f} | benzonase: coverage={cov_t:0.3f} ({nk_t}/{nt_t}), acc_kept={acc_t:0.3f} | "
          f"control: coverage={cov_c:0.3f} ({nk_c}/{nt_c}), acc_kept={acc_c:0.3f}")


print("\nReject option using p_top1")
for thr in [0.10, 0.20, 0.30, 0.40, 0.50]:
    cov_t, acc_t, nk_t, nt_t, _ = reject_metrics(y_sink_treatment, top1_label_t, p1_t, thr)
    cov_c, acc_c, nk_c, nt_c, _ = reject_metrics(y_sink_control,   top1_label_c, p1_c, thr)

    print(f"  thr={thr:0.2f} | benzonase: coverage={cov_t:0.3f} ({nk_t}/{nt_t}), acc_kept={acc_t:0.3f} | "
          f"control: coverage={cov_c:0.3f} ({nk_c}/{nt_c}), acc_kept={acc_c:0.3f}")

print("Evaluated predictions.")

#  9. Figures 
print("Step 9: Generating figures...")

# --- Confusion matrices (normalized by true label) ---
print("Plotting confusion matrices...")

def donor_key(label: str):
    s = str(label).strip()
    m = re.search(r'(\d+)$', s)
    return (0, int(m.group(1))) if m else (1, s)

classes_sorted = sorted(classes, key=donor_key)

cm_treatment = confusion_matrix(
    y_sink_treatment, y_pred_treatment, labels=classes_sorted, normalize="true"
)
cm_control = confusion_matrix(
    y_sink_control, y_pred_control, labels=classes_sorted, normalize="true"
)

plt.figure(figsize=(10, 8))
sns.heatmap(
    cm_treatment, annot=False, cmap="Blues",
    xticklabels=classes_sorted, yticklabels=classes_sorted
)
plt.ylabel("True label")
plt.xlabel("Predicted label")
plt.title("Normalized Confusion Matrix — PVC Benzonase (row-normalized)")
plt.tight_layout()
plt.savefig(f"{OUTPUT}/{MODEL}_pvc_benzonase_confusion_matrix_normalized.pdf", bbox_inches="tight")
plt.close()

plt.figure(figsize=(10, 8))
sns.heatmap(
    cm_control, annot=False, cmap="Blues",
    xticklabels=classes_sorted, yticklabels=classes_sorted
)
plt.ylabel("True label")
plt.xlabel("Predicted label")
plt.title("Normalized Confusion Matrix — PVC Control (row-normalized)")
plt.tight_layout()
plt.savefig(f"{OUTPUT}/{MODEL}_pvc_control_confusion_matrix_normalized.pdf", bbox_inches="tight")
plt.close()

# --- Top feature importances ---
print("Plotting top feature importances...")

top_n = 10
top_features = importance_df.head(top_n).copy()

plt.figure(figsize=(8, 6))
plt.barh(top_features["Feature"], top_features["Importance"])
plt.gca().invert_yaxis()
plt.xlabel("Feature Importance")
plt.title(f"Top {top_n} Features — {MODEL}")
plt.tight_layout()
plt.savefig(f"{OUTPUT}/{MODEL}_top{top_n}_features.pdf", bbox_inches="tight")
plt.close()

# Top-K curve
def topk_curve(y_true, y_proba, classes, k_max=None):
    """
    Returns ks (1..k_max) and top-k hit rates.
    y_true: array-like of true labels (strings)
    y_proba: (n_samples, n_classes) proba aligned to classes
    classes: array-like of class labels (strings) in proba column order
    """
    y_true = np.asarray(y_true, dtype=str)
    classes = np.asarray(classes, dtype=str)

    n_classes = y_proba.shape[1]
    if k_max is None:
        k_max = n_classes
    k_max = min(k_max, n_classes)

    # topk_idx_desc[:, :k] gives indices of top-k predicted classes
    topk_idx_desc = np.argsort(y_proba, axis=1)[:, ::-1]  # descending
    topk_labels_desc = classes[topk_idx_desc]             # (n_samples, n_classes)

    hits = []
    ks = np.arange(1, k_max + 1)
    for k in ks:
        hit_k = np.mean([y_true[i] in topk_labels_desc[i, :k] for i in range(len(y_true))])
        hits.append(float(hit_k))
    return ks, np.array(hits, dtype=float)

k_max = 25  # or 25 for all donors

ks_c, hits_c = topk_curve(y_sink_control,   y_proba_control,   proba_classes, k_max=k_max)
ks_t, hits_t = topk_curve(y_sink_treatment, y_proba_treatment, proba_classes, k_max=k_max)

def save_topk_csv(ks, hits, out_csv, model_name, condition_name):
    df = pd.DataFrame({
        "k": ks,
        "hit_rate": hits,
        "model": model_name,
        "condition": condition_name,
    })
    df.to_csv(out_csv, index=False)

save_topk_csv(ks_c, hits_c, f"{OUTPUT}/{MODEL}_topk_curve_control.csv", f"{MODEL}", "control")
save_topk_csv(ks_t, hits_t, f"{OUTPUT}/{MODEL}_topk_curve_treatment.csv", f"{MODEL}", "benzonase")

plt.figure(figsize=(6, 4))
plt.plot(ks_c, hits_c, marker="o", markersize=3, label="Control")
plt.plot(ks_t, hits_t, marker="o", markersize=3, label="Treatment")
plt.xticks(ks_c, rotation=90)
plt.ylim(0, 1.0)
plt.xlabel("k")
plt.ylabel("Top-k hit rate")
plt.title(f"Top-k curve (source: {MODEL})")
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUTPUT}/{MODEL}_topk_curve_control_vs_treatment.pdf", bbox_inches="tight")
plt.close()

print("Generated figures.")

####################################################################################################################################

# Permutation test (global null model test)
if args.run_permtest:
    print("Step 9: Starting permutation test...")
    logfile.write("Step 9: Starting permutation test...\n")
    # In train mode, best_params is computed after GridSearchCV.
    # In predict mode, best_params must come from the bundle.
    if best_params is None:
        raise ValueError(
            "Permutation test requires best_params. "
            "Run --mode train first (to compute/save best_params), or use predict mode with a bundle that includes best_params."
        )

    # Ensure we reuse the SAME leakage-controlled CV splits as training
    try:
        _ = cv_splits
    except NameError:
        raise ValueError(
            "Permutation test requires cv_splits to be defined (your donor-stratified grouped splits). "
            "Make sure cv_splits is created before this section."
        )

    def permute_labels_global(y, rng):
        """Global permutation of labels (breaks X↔y association)."""
        return rng.permutation(np.asarray(y))

    def cv_balacc_fixed_params(X, y, cv_splits, xgb_params):
        fold_scores = []

        is_df = hasattr(X, "iloc")  # crude but works for pandas

        for train_idx, test_idx in cv_splits:
            if is_df:
                X_tr = X.iloc[train_idx, :].values
                X_te = X.iloc[test_idx, :].values
            else:
                X_tr = X[train_idx]
                X_te = X[test_idx]

            y_tr = np.asarray(y)[train_idx]
            y_te = np.asarray(y)[test_idx]

            le = LabelEncoder()
            y_tr_enc = le.fit_transform(y_tr)
            labels_fold = np.arange(len(le.classes_))

            y_te_mapped = np.array([le.transform([v])[0] if v in le.classes_ else -1 for v in y_te])
            keep = y_te_mapped != -1
            if np.sum(keep) == 0:
                continue

            clf = XGBClassifier(
                objective="multi:softprob",
                eval_metric="mlogloss",
                random_state=42,
                n_jobs=1,
                **xgb_params
            )
            clf.fit(X_tr, y_tr_enc)
            y_pred = clf.predict(X_te[keep])
            fold_scores.append(
                recall_score(
                    y_te_mapped[keep],
                    y_pred,
                    labels=labels_fold,
                    average="macro",
                    zero_division=0
                )
            )

        return float(np.mean(fold_scores)) if fold_scores else np.nan


    # 1) Observed score on real labels
    obs_mean = cv_balacc_fixed_params(X_source, y_source_enc, cv_splits, best_params)

    print("\nPermutation test setup (global label shuffle):")
    print(f"Observed CV balanced_accuracy mean = {obs_mean:.6f}")
    logfile.write("\nPermutation test setup (global label shuffle):\n")
    logfile.write(f"Observed CV balanced_accuracy mean = {obs_mean:.6f}\n")

    # 2) Permutations
    N_PERM = 999
    rng = np.random.default_rng(2026)
    perm_means = np.empty(N_PERM, dtype=float)

    for i in range(N_PERM):
        y_perm = permute_labels_global(y_source_enc, rng)
        perm_means[i] = cv_balacc_fixed_params(X_source, y_perm, cv_splits, best_params)
        if (i + 1) % 50 == 0:
            print(f"  completed {i+1}/{N_PERM} permutations")

    # 3) One-sided p-value: P(null >= observed) with +1 smoothing
    p_value = (1 + np.sum(perm_means >= obs_mean)) / (N_PERM + 1)

    print("\nPermutation test results (global label shuffle):")
    print(f"Observed mean bal_acc = {obs_mean:.6f}")
    print(f"Null mean bal_acc     = {perm_means.mean():.6f}")
    print(f"Null max bal_acc      = {perm_means.max():.6f}")
    print(f"p-value (>= obs)      = {p_value:.6g}")

    logfile.write("\nPermutation test results (global label shuffle):\n")
    logfile.write(f"Observed mean bal_acc = {obs_mean:.6f}\n")
    logfile.write(f"Null mean bal_acc     = {perm_means.mean():.6f}\n")
    logfile.write(f"Null max bal_acc      = {perm_means.max():.6f}\n")
    logfile.write(f"p-value (>= obs)      = {p_value:.6g}\n")

    # 4) Save null distribution
    perm_df = pd.DataFrame({"perm_mean_balacc": perm_means})
    perm_out = f"{OUTPUT}/{MODEL}_permtest_balacc.csv"
    perm_df.to_csv(perm_out, index=False)
    print(f"Saved permutation null distribution to {perm_out}")
    logfile.write(f"Saved permutation null distribution to {perm_out}\n")

    # 5) Plot histogram
    plt.figure(figsize=(7, 5))
    plt.hist(perm_means, bins=30)
    plt.axvline(obs_mean, linewidth=2)
    plt.xlabel("Mean CV balanced accuracy (permuted labels)")
    plt.ylabel("Count")
    plt.title("Permutation test null distribution")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT}/{MODEL}_permtest_balacc_hist.pdf", bbox_inches="tight")
    plt.close()


logfile.close()