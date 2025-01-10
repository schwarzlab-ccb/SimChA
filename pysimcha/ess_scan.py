#%%
import numpy as np
from numba import njit
import pandas as pd
import matplotlib.pyplot as plt
import os
from os.path import join
import seaborn as sns
from collections import defaultdict
# %%
p_hemi = 0.184364
p_null = 0.000308
p_rest = 1.0 - p_hemi - p_null
pcawg_freq = np.array([p_hemi, p_null, p_rest])
print(pcawg_freq)
def err_func(sim_freq):
    return -np.sum(pcawg_freq * np.log(sim_freq[:len(pcawg_freq)]))

# %%
sim_dir = "../out/iscb_results/results_all_events_noWGD"
subdirs = [d for d in os.listdir(sim_dir) if os.path.isdir(join(sim_dir, d))]
sim_res = []
for i, subdir in enumerate(subdirs):
    print(f"Folders processed: {i}/{len(subdirs)}")
    ess = float(subdir)
    clones_df = pd.read_csv(join(sim_dir, subdir, "clones.tsv"), sep="\t")
    hemi_samples = np.mean(clones_df.hemizygosity.values)
    null_samples = np.mean(clones_df.nullizygosity.values)
    rest_samples = 1.0 - hemi_samples - null_samples
    err = err_func(np.array([hemi_samples, null_samples, rest_samples]))
    sim_res.append([ess, err, hemi_samples, null_samples, rest_samples])
    #ploidy = clones_df.ploidy.values
    #mean_ploidy = np.mean(ploidy)
    #var_ploidy = np.var(ploidy)
    #err = err_func(mean_ploidy)
    #abs_err = abs(err)
    #sim_res.append([ess, err, abs_err, mean_ploidy, var_ploidy])

# %%
df = pd.DataFrame(sim_res, columns=["ess", "error", "hemi", "null", "rest"])
df.to_csv("ess_scan.tsv", sep="\t", index=False)
#df_2 = pd.DataFrame(sim_res_2, columns=["ess", "error", "abs_err", "mean_ploidy", "var_ploidy"])
plt.figure(figsize=(10, 8))
plt.scatter(df["ess"], df["error"], alpha=1)
plt.xlabel(r"Essentiality parameter - $\gamma$")
plt.ylabel("Log-Likelihood")
plt.title(r"Effect of Essentiality on WGD- ploidies")
plt.savefig(join("../img","diploid_ess_scan.png"), dpi=450, bbox_inches="tight")
#%%
print(f"PCAWG:\t hemi: {pcawg_freq[0]:6f}; nulli: {pcawg_freq[1]:6f}; rest: {pcawg_freq[2]:6f}")
best_error = df.query("error == error.min()")
print(f"Best:\t hemi: {best_error.hemi.values[0]:6f}; nulli: {best_error.null.values[0]:6f}; rest: {best_error.rest.values[0]:6f}")
worst_error = df.query("error == error.max()")
print(f"Worst:\t hemi: {worst_error.hemi.values[0]:6f}; nulli: {worst_error.null.values[0]:6f}; rest: {worst_error.rest.values[0]:6f}")
