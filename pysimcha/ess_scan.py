#%%
import numpy as np
from numba import njit
import pandas as pd
import matplotlib.pyplot as plt
import os
from os.path import join
# %%
data_dir = "../out/pcawg"
pcawg_clones = pd.read_csv(join(data_dir, "clones.tsv"), sep="\t")
nonwgd_samples = pd.read_csv(join("../data", "pcawg_nonWGD_event_counts.tsv"), sep="\t")
# %%
nonwgd_ids = nonwgd_samples.sample_id.values
nonwgd_ploidies = pcawg_clones[pcawg_clones["sample_id"].isin(nonwgd_ids)].ploidy.values
print(np.mean(nonwgd_ploidies))
# %%
def err_func(mean_ploidy):
    return (mean_ploidy - 2)

# %%
sim_dir = "../out/nf/results_dyn_mut_ess_scan"
subdirs = [d for d in os.listdir(sim_dir) if os.path.isdir(join(sim_dir, d))]
sim_res = []
selected_dirs = [0.75, 0.8, 0.85, 0.88, 0.936]
selected_res = {d : [] for d in selected_dirs}
all_ploidies = {d : [] for d in selected_dirs}
for i, subdir in enumerate(subdirs):
    print(f"Folders processed: {i}/{len(subdirs)}")
    ess = float(subdir)
    clones_df = pd.read_csv(join(sim_dir, subdir, "clones.tsv"), sep="\t")
    events_df = pd.read_csv(join(sim_dir, subdir, "events.tsv"), sep="\t")
    wgd_samples = events_df.query("event_type == 'WholeGenomeDoubling'").sample_id.unique()
    ploidy = clones_df[~clones_df["sample_id"].isin(wgd_samples)].ploidy.values
    mean_ploidy = np.mean(ploidy)
    var_ploidy = np.var(ploidy)
    err = err_func(mean_ploidy)
    abs_err = abs(err)
    sim_res.append([ess, err, abs_err, mean_ploidy, var_ploidy])
    if ess in selected_dirs:
        selected_res[ess] = ploidy
        all_ploidies[ess] = clones_df.ploidy.values
#%%
sim_dir = "../out/nf/results_ess_scan_mut_0.5_zoomed"
subdirs = [d for d in os.listdir(sim_dir) if os.path.isdir(join(sim_dir, d))]
sim_res_2 = []
for i, subdir in enumerate(subdirs):
    print(f"Folders processed: {i}/{len(subdirs)}")
    ess = float(subdir)
    clones_df = pd.read_csv(join(sim_dir, subdir, "clones.tsv"), sep="\t")
    events_df = pd.read_csv(join(sim_dir, subdir, "events.tsv"), sep="\t")
    wgd_samples = events_df.query("event_type == 'WholeGenomeDoubling'").sample_id.unique()
    ploidy = clones_df[~clones_df["sample_id"].isin(wgd_samples)].ploidy.values
    mean_ploidy = np.mean(ploidy)
    var_ploidy = np.var(ploidy)
    err = err_func(mean_ploidy)
    abs_err = abs(err)
    sim_res_2.append([ess, err, abs_err, mean_ploidy, var_ploidy])
# %%
df = pd.DataFrame(sim_res, columns=["ess", "error", "abs_err", "mean_ploidy", "var_ploidy"])
df.to_csv("ess_scan.tsv", sep="\t", index=False)
#df_2 = pd.DataFrame(sim_res_2, columns=["ess", "error", "abs_err", "mean_ploidy", "var_ploidy"])
plt.figure(figsize=(10, 8))
plt.scatter(df["ess"], df["error"], alpha=0.8, label="Mut Rate = 0.25")
for i, ess in enumerate(selected_dirs):
    val = df.query(f"ess == {ess}").error.iloc[0]
    plt.scatter(ess, val, color=f"C{i+2}", marker="x", s=150, lw=4, label=f"Ess: {ess}")

#plt.scatter(df_2["ess"], df_2["error"], alpha=0.8, label="Mut Rate = 0.5")
plt.xlabel("Essentiality")
plt.ylabel("Ploidy difference from diploid")
plt.title("Essentiality effect on ploidy with dynamic mutation rate")
plt.legend()
plt.savefig("dyn_mut_ess_scan", dpi=300, bbox_inches="tight")
#%%
dir = "../out/nf"
clones_normal = pd.read_csv(join(dir, "results_double_ess_scan_0.93_more_loss", "clones.tsv"), sep="\t")
clones_dyn = pd.read_csv(join(dir, "results_triple_ess_scan_0.93_more_loss", "clones.tsv"), sep="\t")
ploidy_normal = clones_normal.ploidy.values
ploidy_dyn = clones_dyn.ploidy.values

plt.figure(figsize=(10, 8))
bins = np.linspace(0, 8, 100)
plt.hist(pcawg_clones.ploidy, bins=bins, histtype="step", label="PCAWG samples", color="black", density=True)
plt.hist(ploidy_normal, bins=bins, histtype="step", label="Ploidy-indep.", color="red", density=True)
plt.hist(ploidy_dyn, bins=bins, histtype="step", label="Ploidy-dep.", color="blue", density=True)
plt.xlabel("Ploidy")
plt.ylabel("Density")
plt.title("Ploidy distribution for Essentiality = 0.9")
plt.legend()
plt.savefig("dyn_mut_ploidy.png", dpi=300, bbox_inches="tight")
# %%
