import os
import glob
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

from skimage import io, measure, morphology
from scipy import ndimage
from scipy.spatial.distance import cdist

# %% ================= USER INPUT =================

root_path = r'F:\Lab\MouseMusclePOGLUT1\results\detectionMetrics'

gt_folder = os.path.join(root_path, 'imagesTest', 'GT emi_vicky combined')

models = {
    'modelV1': os.path.join(root_path, 'imagesTest', 'modelV1'),
    'modelV2': os.path.join(root_path, 'imagesTest', 'modelV2')
}

dist_thresh = 10
img_ext = '*.tif'

microns_per_pixel = 0.5055443
min_area_um2 = 10
min_area_pixels = math.ceil(min_area_um2 / (microns_per_pixel**2))

thresholds_pct = np.arange(0.30, 1.01, 0.02)
thresholds = np.unique((thresholds_pct * 255).astype(np.uint8))
thresholds_pct = thresholds / 255

# %% ================= OUTPUT FOLDER =================

date_str = datetime.now().strftime('%Y%m%d')
out_dir = os.path.join(root_path, f"threshold_optimization_{date_str}")
os.makedirs(out_dir, exist_ok=True)

# %% ================= LOAD GT =================

gt_files = glob.glob(os.path.join(gt_folder, img_ext))
eps = np.finfo(float).eps

gt_data = {}
for img_path in gt_files:
    file_name = os.path.basename(img_path)
    gt_data[file_name] = io.imread(img_path)

print(f"Loaded {len(gt_data)} GT images")

# %% ================= FUNCTION =================

def evaluate_model(pred_folder, model_name):

    print(f"\n=== Evaluating {model_name} ===")

    data = []
    for file_name, gt in gt_data.items():
        pr = io.imread(os.path.join(pred_folder, file_name))

        if pr.dtype != np.uint8:
            BW = (pr * 255).astype(np.uint8)
        else:
            BW = pr

        data.append((file_name, gt, BW))

    results = []

    for thr, thr_pct in zip(thresholds, thresholds_pct):

        Precisions, Recalls, F1s = [], [], []

        TP_total, FP_total, FN_total = 0, 0, 0  # 👈 micro accumulators

        for file_name, gt, BW in data:

            gt_bin = gt > 0
            pr_bin = BW > thr

            pr_bin = ndimage.binary_fill_holes(pr_bin)
            pr_bin = morphology.binary_closing(pr_bin, morphology.disk(1))
            pr_bin = morphology.remove_small_objects(pr_bin, min_size=min_area_pixels)

            gt_L = measure.label(gt_bin, connectivity=2)
            pr_L = measure.label(pr_bin, connectivity=2)

            gt_props = measure.regionprops(gt_L)
            pr_props = measure.regionprops(pr_L)

            n_GT = len(gt_props)
            n_PR = len(pr_props)

            if n_GT == 0 and n_PR == 0:
                Precisions.append(1); Recalls.append(1); F1s.append(1)
                continue
            if n_GT == 0:
                Precisions.append(0); Recalls.append(1); F1s.append(0)
                FP_total += n_PR
                continue
            if n_PR == 0:
                Precisions.append(1); Recalls.append(0); F1s.append(0)
                FN_total += n_GT
                continue

            gt_C = np.array([p.centroid for p in gt_props])
            pr_C = np.array([p.centroid for p in pr_props])

            D = cdist(gt_C, pr_C)

            matched_GT = np.zeros(n_GT, dtype=bool)
            matched_PR = np.zeros(n_PR, dtype=bool)

            TP = 0
            while True:
                min_D = np.min(D)
                if min_D > dist_thresh:
                    break

                g, p = np.unravel_index(np.argmin(D), D.shape)
                TP += 1
                matched_GT[g] = True
                matched_PR[p] = True

                D[g, :] = np.inf
                D[:, p] = np.inf

            FN = np.sum(~matched_GT)
            FP = np.sum(~matched_PR)

            TP_total += TP
            FP_total += FP
            FN_total += FN

            p_val = TP / (TP + FP + eps)
            r_val = TP / (TP + FN + eps)
            f1_val = 2 * p_val * r_val / (p_val + r_val + eps)

            Precisions.append(p_val)
            Recalls.append(r_val)
            F1s.append(f1_val)

        # Macro
        precision_macro = np.mean(Precisions)
        recall_macro = np.mean(Recalls)
        f1_macro = np.mean(F1s)

        # Micro
        precision_micro = TP_total / (TP_total + FP_total + eps)
        recall_micro = TP_total / (TP_total + FN_total + eps)
        f1_micro = 2 * precision_micro * recall_micro / (precision_micro + recall_micro + eps)

        results.append({
            'Model': model_name,
            'Threshold': thr,
            'Threshold_%': thr_pct,
            'Precision_macro': precision_macro,
            'Recall_macro': recall_macro,
            'F1_macro': f1_macro,
            'Precision_micro': precision_micro,
            'Recall_micro': recall_micro,
            'F1_micro': f1_micro
        })

    df = pd.DataFrame(results)

    out_excel = os.path.join(out_dir, f"{model_name}_metrics_macro_micro.xlsx")
    df.to_excel(out_excel, index=False)

    print(f"Saved {model_name} results to: {out_excel}")

    return df

# %% ================= RUN =================

all_results = []

for model_name, pred_folder in models.items():
    df_model = evaluate_model(pred_folder, model_name)
    all_results.append(df_model)

df_all = pd.concat(all_results, ignore_index=True)

combined_excel = os.path.join(out_dir, "combined_models_macro_micro.xlsx")
df_all.to_excel(combined_excel, index=False)

print(f"\nSaved combined results to: {combined_excel}")

# %% ================= BEST =================

for model_name in models.keys():
    df_m = df_all[df_all['Model'] == model_name]

    best_macro = df_m.loc[df_m['F1_macro'].idxmax()]
    best_micro = df_m.loc[df_m['F1_micro'].idxmax()]

    print(f"\n{model_name}")
    print("Best MACRO F1:")
    print(best_macro[['Threshold_%', 'F1_macro']])

    print("Best MICRO F1:")
    print(best_micro[['Threshold_%', 'F1_micro']])

# %% ================= PLOT =================

plt.figure()

for model_name in models.keys():
    df_m = df_all[df_all['Model'] == model_name]

    plt.plot(df_m['Threshold_%']*100, df_m['F1_macro'],
             label=f'{model_name} Macro F1')

    plt.plot(df_m['Threshold_%']*100, df_m['F1_micro'],
             linestyle='--', label=f'{model_name} Micro F1')

plt.xlabel('Threshold (% of max intensity)')
plt.ylabel('F1 Score')
plt.title('Macro vs Micro F1 Comparison')
plt.legend()
plt.grid(True)

plot_path = os.path.join(out_dir, "macro_vs_micro_comparison.png")
plt.savefig(plot_path, dpi=300)
plt.show()

print(f"Saved plot to: {plot_path}")