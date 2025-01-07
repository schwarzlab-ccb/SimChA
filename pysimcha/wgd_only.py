# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from os.path import join

# %%
def get_synthetic_counts(dir):
    clones_df = pd.read_csv(join(dir, "clones.tsv"), sep="\t")
    n = len(clones_df)
    events_df = pd.read_csv(join(dir, "events.tsv"), sep="\t")
    wgd_events = events_df[events_df["event_type"] == "WholeGenomeDoubling"]
    wgd_counts = wgd_events["sample_id"].value_counts()
    zero_count = n - wgd_events.sample_id.nunique()
    counts = np.ones(4)
    for count in wgd_counts:
        if count >= len(counts):
            counts[-1] += 1
        else:
            counts[count] += 1
    counts[0] += zero_count
    counts = np.array(counts)
    return counts / np.sum(counts)

pcawg_counts = np.array([1960, 818, 12, 0])
pcawg_counts += 1
pcawg_freq = pcawg_counts / np.sum(pcawg_counts)
def log_like(sim_freq):
    return -np.sum(pcawg_freq * np.log(sim_freq))
# %%
out_dir = "../out/iscb_results/results_wgd_and_alpha"
subdirs = [d for d in os.listdir(out_dir) if os.path.isdir(join(out_dir, d))]
simcha_res = []
for subdir in subdirs:
    counts = get_synthetic_counts(join(out_dir, subdir))
    params = subdir.split("_")
    alpha = float(params[1])
    p_wgd = float(params[0])
    simcha_res.append([p_wgd, alpha, log_like(counts)])
# %%
df = pd.DataFrame(simcha_res, columns=["p_wgd", "alpha", "error"])
pivot_df = df.pivot_table(index="p_wgd", columns="alpha", values="error")

# Find the minimum error value and its corresponding p_wgd and alpha
min_error = pivot_df.min().min()
min_error_location = np.unravel_index(np.argmin(pivot_df.values), pivot_df.shape)
best_p_wgd = pivot_df.index[min_error_location[0]]
best_alpha = pivot_df.columns[min_error_location[1]]

# Plot the error surface
plt.figure(figsize=(10, 8))
contour = plt.contourf(pivot_df.columns, pivot_df.index, pivot_df.values, levels=100, cmap='viridis', vmax=3)
plt.colorbar(contour)
label = r'Best score: $\alpha$ - ' +f"{best_alpha:.3f}, " +r"$p_\text{WGD}$ - " + f"{best_p_wgd:.3f}"
plt.plot(best_alpha, best_p_wgd, 'r*', markersize=12, label=label)
plt.xlabel(r'Stress parameter - $\alpha$')
plt.ylabel(r'$p_\text{WGD}$')
plt.title(f'Stress vs Whole Genome Doubling Error Surface')
plt.legend()
plt.savefig(join("../img","wgd_only_dyn_mut_rate_4.png"), dpi=450, bbox_inches="tight")
plt.show()
# %%
out_dir = "../out/iscb_results/results_wgd_and_alpha_no_dyn"
subdirs = [d for d in os.listdir(out_dir) if os.path.isdir(join(out_dir, d))]
simcha_res = []
for subdir in subdirs:
    counts = get_synthetic_counts(join(out_dir, subdir))
    params = subdir.split("_")
    alpha = float(params[1])
    p_wgd = float(params[0])
    simcha_res.append([p_wgd, alpha, log_like(counts)])
# %%
df = pd.DataFrame(simcha_res, columns=["p_wgd", "alpha", "error"])
pivot_df = df.pivot_table(index="p_wgd", columns="alpha", values="error")

# Find the minimum error value and its corresponding p_wgd and alpha
min_error = pivot_df.min().min()
min_error_location = np.unravel_index(np.argmin(pivot_df.values), pivot_df.shape)
best_p_wgd = pivot_df.index[min_error_location[0]]
best_alpha = pivot_df.columns[min_error_location[1]]

# Plot the error surface
plt.figure(figsize=(10, 8))
contour = plt.contourf(pivot_df.columns, pivot_df.index, pivot_df.values, levels=100, cmap='viridis', vmax=3)
plt.colorbar(contour)
label = r'Best score: $\alpha$ - ' +f"{best_alpha:.3f}, " +r"$p_\text{WGD}$ - " + f"{best_p_wgd:.3f}"
plt.plot(best_alpha, best_p_wgd, 'r*', markersize=12, label=label)
plt.xlabel(r'Stress parameter - $\alpha$')
plt.ylabel(r'$p_\text{WGD}$')
plt.title(f'Stress vs Whole Genome Doubling Error Surface (Without Dynamic Mutation Rate)')
plt.legend()
plt.savefig(join("../img","wgd_only_no_dyn_mut_rate.png"), dpi=450, bbox_inches="tight")
plt.show()
# %%
