# %%
import numpy as np
from numba import njit
import pandas as pd
import matplotlib.pyplot as plt
import os
from os.path import join
# %%
data_dir = "../out/pcawg"
pcawg_clones = pd.read_csv(join(data_dir, "clones.tsv"), sep="\t")
wgd_samples = pd.read_csv(join("../data", "pcawg_WGD_event_counts.tsv"), sep="\t")
# %%
# Get all the ploidies of the WGD-positive samples
wgd_ids = wgd_samples.sample_id.values
wgd_ploidies = pcawg_clones[pcawg_clones["sample_id"].isin(wgd_ids)].ploidy.values
# Remove any samples with ploidy > 5
wgd_ploidies = wgd_ploidies[wgd_ploidies <= 5]
mean_pcawg_ploidy = np.mean(wgd_ploidies)
var_pcawg_ploidies = np.var(wgd_ploidies)
print(mean_pcawg_ploidy)

#%%
def err_func(mean_ploidy):
    return (mean_ploidy - mean_pcawg_ploidy)**2
# %%
sim_dir = "../out/nf/results_loss_theta_20"
subdirs = [d for d in os.listdir(sim_dir) if os.path.isdir(join(sim_dir, d))]
sim_res = []
for subdir in subdirs:
    params = subdir.split("_")
    r_loss = float(params[0])
    alpha = float(params[1])
    clones_df = pd.read_csv(join(sim_dir, subdir, "clones.tsv"), sep="\t")
    ploidy = clones_df.ploidy.values
    mean_ploidy = np.mean(ploidy)
    var_ploidy = np.var(ploidy)
    err = err_func(mean_ploidy)

    sim_res.append([r_loss, alpha, err, mean_ploidy, var_ploidy])
#%%
df = pd.DataFrame(sim_res, columns=["r_loss", "alpha", "error", "mean_ploidy", "var_ploidy"])
pivot_df = df.pivot_table(index="r_loss", columns="alpha", values="error")
# Find the minimum error value and its corresponding p_wgd and alpha
min_error = pivot_df.min().min()
min_error_location = np.unravel_index(np.argmin(pivot_df.values), pivot_df.shape)
best_r_loss = pivot_df.index[min_error_location[0]]
best_alpha = pivot_df.columns[min_error_location[1]]

# Plot the error surface
plt.figure(figsize=(10, 8))
contour = plt.contourf(pivot_df.columns, pivot_df.index, pivot_df.values, levels=100, cmap='viridis')
plt.colorbar(contour)
plt.plot(best_alpha, best_r_loss, 'r*', markersize=12, label='Best score')
plt.xlabel('Alpha')
plt.ylabel('r_loss')
plt.title(f'SimChA Error Surface - best r_loss: {best_r_loss}, best alpha: {best_alpha}')
plt.legend()
plt.savefig("simcha_error_surface_losses.png")
plt.show()
# %%
df.to_csv("simcha_loss_scan_results.tsv", sep="\t", index=False)
# %%
