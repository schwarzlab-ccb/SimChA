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
data_dir = "../out/pcawg"
pcawg_clones = pd.read_csv(join(data_dir, "clones.tsv"), sep="\t")
nonwgd_samples = pd.read_csv(join("../data", "pcawg_nonWGD_event_counts.tsv"), sep="\t")
# %%
nonwgd_ids = nonwgd_samples.sample_id.values
nonwgd_ploidies = pcawg_clones[pcawg_clones["sample_id"].isin(nonwgd_ids)].ploidy.values
mean_pcawg_ploidy = np.mean(nonwgd_ploidies)
# %%
p_hemi = 0.184364
p_null = 0.000308
p_rest = 1.0 - p_hemi - p_null
pcawg_freq = np.array([p_hemi, p_null, p_rest])
print(pcawg_freq)
def err_func(sim_freq):
    return -np.sum(pcawg_freq * np.log(sim_freq[:len(pcawg_freq)]))

# %%
sim_dir = "../out/iscb_results/results_diploid_ess_scan_extensive"
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
best_error = df.query("error == error.min()")
print(f"hemi: {best_error.hemi.values[0]:6f}; nulli: {best_error.null.values[0]:6f}; rest: {best_error.rest.values[0]:6f}")
worst_error = df.query("error == error.max()")
print(f"hemi: {worst_error.hemi.values[0]:6f}; nulli: {worst_error.null.values[0]:6f}; rest: {worst_error.rest.values[0]:6f}")

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
dir = "../out/nf"
ls = ["-", "-", "-", "-"]

plt.figure(figsize=(10, 8))
bins = np.linspace(1.4, 5, 50)

plt.hist(pcawg_clones.ploidy, bins=bins, histtype="step", label="PCAWG samples", lw=1.5, color="black", density=True)

"""plt.hist(pd.read_csv(join(dir, "results_wgd_separated_mut_1_0.9", "clones.tsv"), sep="\t").ploidy.values, 
         bins=bins, histtype="step", label=r"$\mu$=1, Haploinsuf.", ls="-.", density=True)
plt.hist(pd.read_csv(join(dir, "results_wgd_separated_mut_1_haplosuf_0.9", "clones.tsv"), sep="\t").ploidy.values, 
         bins=bins, histtype="step", label=r"$\mu$=1", ls="--", density=True)"""

plt.hist(pd.read_csv(join(dir, "results_wgd_separated_0.9", "clones.tsv"), sep="\t").ploidy.values, 
         bins=bins, histtype="step", label=r"$\mu$=0.5, Haploinsuf.", ls="-.", lw=1.5, density=True)
plt.hist(pd.read_csv(join(dir, "results_wgd_separated_haplosuf_0.9", "clones.tsv"), sep="\t").ploidy.values, 
         bins=bins, histtype="step", label=r"$\mu$=0.5", ls="--", lw=1.5, density=True)

plt.xlabel("Ploidy")
plt.ylabel("Density")
plt.title("Ploidy distribution for mut_rate = 0.5, fit params {0.1, 0, 0.9}, Separated WGD-status events")
plt.legend()
plt.savefig("dyn_mut_ploidy_mut_rate_0.5.png", dpi=300, bbox_inches="tight")
# %%
def get_wgd_status(events_file, suffix):
    events = pd.read_csv(events_file, sep="\t")
    # Filter for WholeGenomeDoubling events
    wgd_events = events[events["event_type"] == "WholeGenomeDoubling"]
    # Create a dictionary with sample_id as keys and 1 if WGD event exists, else 0
    wgd_status = {f"{sample_id}_{suffix}": 1 for sample_id in wgd_events["sample_id"].unique()}
    # Add entries for sample_ids without WGD events
    all_sample_ids = events["sample_id"].unique()
    for sample_id in all_sample_ids:
        id = f"{sample_id}_{suffix}"
        if id not in wgd_status:
            wgd_status[id] = 0
    return wgd_status
def plot_average_inter_event_time(events_file, suffix, ax=None):
    if ax is None:
        ax = plt.gca()
    events = pd.read_csv(events_file, sep="\t")
    # Calculate inter-event times
    events["time_diff"] = events.groupby("sample_id")["time"].diff()
    # Filter out the first event (NaN time_diff)
    inter_event_times = events.dropna(subset=["time_diff"])
    # Calculate average inter-event time for each sample_id
    avg_inter_event_times = inter_event_times.groupby("sample_id")["time_diff"].mean()
    # Get WGD status
    wgd_status = get_wgd_status(events_file, suffix)
    # Separate WGD and non-WGD samples
    non_wgd_samples = [avg_inter_event_times[f"{sample_id}_{suffix}"] for sample_id in avg_inter_event_times.index if wgd_status[f"{sample_id}_{suffix}"] == 0]
    wgd_samples = [avg_inter_event_times[f"{sample_id}_{suffix}"] for sample_id in avg_inter_event_times.index if wgd_status[f"{sample_id}_{suffix}"] == 1]
    # Plot using seaborn.violinplot
    sns.violinplot(data=[non_wgd_samples, wgd_samples], ax=ax)
# %%
fig, ax = plt.subplots(1, 1, figsize=(10,10))
plot_average_inter_event_time(join("../out/nf", "results_wgd_separated_0.9", "events.tsv"), "a", ax=ax)

ax.set_title("Continuous time model with ploidy-dependent mutation rate, Fitness & Dynamic Rate")
ax.set_ylabel("Mean Inter-event Time")
ax.set_xlabel("WGD status")
fig.savefig("../img/inter_event_time_wgd_separated.png", dpi=300, bbox_inches="tight")
# %%
dir = "../out/wgd_test_0.5_fit_haplosuf_sqr_stress"
plt.figure(figsize=(10, 8))
bins = np.linspace(1.2, 8, 50)

plt.hist(pcawg_clones.ploidy, bins=bins, histtype="step", label="PCAWG samples", lw=1.5, color="black", density=True)

"""plt.hist(pd.read_csv(join(dir, "results_wgd_separated_mut_1_0.9", "clones.tsv"), sep="\t").ploidy.values, 
         bins=bins, histtype="step", label=r"$\mu$=1, Haploinsuf.", ls="-.", density=True)
plt.hist(pd.read_csv(join(dir, "results_wgd_separated_mut_1_haplosuf_0.9", "clones.tsv"), sep="\t").ploidy.values, 
         bins=bins, histtype="step", label=r"$\mu$=1", ls="--", density=True)"""
clones = pd.read_csv(join(dir, "clones.tsv"), sep="\t")
# Find all the clones with no events
plt.hist(clones.ploidy.values, 
         bins=bins, histtype="step", label=r"$ \mu $=0.1", ls="-.", lw=1.5, density=True)

clones = pd.read_csv(join("../out/wgd_test_0.5_no_fit_haplosuf_sqr_stress", "clones.tsv"), sep="\t")
# Find all the clones with no events
plt.hist(clones.ploidy.values, 
         bins=bins, histtype="step", label=r"$ \mu $=0.1, No Fitness", ls="-.", lw=1.5, density=True)

plt.xlabel("Ploidy")
plt.ylabel("Density")
plt.title("Ploidy distribution with large-scale events simulated only")
plt.legend()
plt.savefig("fit_no_fit_comp_zoomed_out.png", dpi=300, bbox_inches="tight")
# %%

plt.figure(figsize=(10, 8))
bins = np.linspace(1.25, 8, 50)

plt.hist(pcawg_clones.ploidy, bins=bins, histtype="step", label="PCAWG samples", lw=1.5, color="black", density=True)

"""plt.hist(pd.read_csv(join(dir, "results_wgd_separated_mut_1_0.9", "clones.tsv"), sep="\t").ploidy.values, 
         bins=bins, histtype="step", label=r"$\mu$=1, Haploinsuf.", ls="-.", density=True)
plt.hist(pd.read_csv(join(dir, "results_wgd_separated_mut_1_haplosuf_0.9", "clones.tsv"), sep="\t").ploidy.values, 
         bins=bins, histtype="step", label=r"$\mu$=1", ls="--", density=True)"""
dir = "../out/optimized_losses"
clones = pd.read_csv(join(dir, "clones.tsv"), sep="\t")
# Find all the clones with no events
plt.hist(clones.ploidy.values, 
         bins=bins, histtype="step", label="Optimized with Haploinsuf", ls="-.", lw=1.5, density=True)

dir = "../out/optimized_losses_haplosuf"
clones = pd.read_csv(join(dir, "clones.tsv"), sep="\t")
# Find all the clones with no events
plt.hist(clones.ploidy.values, 
         bins=bins, histtype="step", label="Optimized with Haplosuf", ls="-.", lw=1.5, density=True)
plt.xlabel("Ploidy")
plt.ylabel("Density")
plt.title("Ploidy distribution of Optimized Loss, Alpha, and Delta scan, w/ and w/o haplosufficiency")
plt.legend()
plt.savefig("optimized_losses.png", dpi=300, bbox_inches="tight")
# %%
import matplotlib.pyplot as plt
from cns.utils.selection import cns_head
from cns.display.plot import fig_lines, fig_CN_heatmap, fig_dots
from cns.process.binning import group_bins
from cns.process.binning import add_cns_loc, sum_cns
from cns.utils.selection import cns_head
from cns.utils.files import load_cns, load_samples
from cns.utils.assemblies import get_assembly
from os.path import join
from cns.data_utils import out_path
hg19 = get_assembly("hg19")


# %%
def get_binned_data(path):
    data = load_cns(path)
    data = sum_cns(add_cns_loc(data))
    return data

binned_pcawg_cnps = get_binned_data(join(out_path, "PCAWG_bin_3MB.tsv"))
binned_tcga_cnps = get_binned_data(join(out_path, "TCGA_hg19_bin_3MB.tsv"))
fig_lines([group_bins(binned_pcawg_cnps), group_bins(binned_tcga_cnps)], ["PCAWG", "TCGA - hg19"], column="total_cn")
plt.title("Copy-Number Profiles of PCAWG and TCGA - hg19")
plt.savefig("pan_cancer_cn_profiles.png", dpi=300, bbox_inches='tight')
# %%
