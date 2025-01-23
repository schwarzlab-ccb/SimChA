# %%
import numpy as np
from scipy.stats import wasserstein_distance
import pandas as pd
import matplotlib.pyplot as plt
import cns
import os
from os.path import join as pjoin
import subprocess
import sys

wgd_status = [0, 1]

base_dir = "/projects/ag-schwarzr/project-simcha/simcha"
tetraploid_base_dir = pjoin(base_dir, "results_ISMB_tsg_tetraploid_and_delta_scan")
diploid_base_dir = pjoin(base_dir, "results_ISMB_tsg_diploid_and_delta_scan")

subdir_dict = {}
for status, curr_dir in enumerate([diploid_base_dir, tetraploid_base_dir]):
    subdirs = [d for d in os.listdir(curr_dir) if os.path.isdir(pjoin(curr_dir, d))]
    subdir_dict[status] = subdirs

group_obs_dict = {}

for status in wgd_status:
    group_obs_dict[status] = pd.read_csv(pjoin(base_dir, "data",f"grouped_obs_cnps_wgd_{status}.tsv"), sep="\t")
# %%
def get_cnps(subdir, status):
    curr_dir = tetraploid_base_dir if status == 1 else diploid_base_dir
    df = cns.load_cns(pjoin(curr_dir, subdir, "cns_3MB.tsv"))
    df.query("chrom != 'chrX' and chrom != 'chrY'", inplace=True)
    cnps = cns.group_samples(df)
    cnps["total_cn"] = cnps["hap1_cn"] + cnps["hap2_cn"]
    return cnps

# %%
df_1 = pd.read_csv(pjoin(base_dir, "out", "tsgog_tet_results_mean.tsv"), sep="\t")
df_0 = pd.read_csv(pjoin(base_dir, "out", "tsgog_dip_results_mean.tsv"), sep="\t")
results_df = pd.concat([df_1, df_0])#pd.read_csv("tsg_og_optimization_results.tsv", sep="\t")
wgd_1_df = results_df.query("status == 1")
wgd_0_df = results_df.query("status == 0")

best_wgd_1_result = wgd_1_df.query("wd == wd.min()")
beta_1  = best_wgd_1_result.beta.values[0]
delta_1 = best_wgd_1_result.delta.values[0]
cns_1 = get_cnps(f"{int(beta_1)}_{delta_1}", 1)
#cns_1 = groups_df[1][f"{beta_1}_{delta_1}"]
cns_1["sample_id"] = "Tet. - " + r"$\beta$: " + f"{beta_1:.2f}, " + r"$\delta$: "+ f"{delta_1:.2f}"
best_wgd_0_result = wgd_0_df.query("wd == wd.min()")
beta_0  = best_wgd_0_result.beta.values[0]
delta_0 = best_wgd_0_result.delta.values[0]
cns_0 = get_cnps(f"{int(beta_0)}_{delta_0}", 0)
#cns_0 = groups_df[0][f"{beta_0}_{delta_0}"]
cns_0["sample_id"] = "Dip. - " + r"$\beta$: " + f"{beta_0:.2f}, " + r"$\delta$: "+ f"{delta_0:.2f}"

d = {"Obs WGD+": group_obs_dict[1], "Obs WGD-": group_obs_dict[0], "Tet. Syn.": cns_1, "Diploid Syn.": cns_0}
combined_df = pd.concat(d.values())
cns.fig_lines(combined_df, cn_columns="total_cn")
plt.savefig(pjoin(base_dir, "img/cnps_tsg_best.png"), dpi=300, bbox_inches="tight")
plt.savefig(pjoin(base_dir, "img/cnps_tsg_best.pdf"), bbox_inches="tight")

# %%
pivot_df = wgd_1_df.pivot_table(index="delta", columns="beta", values="wd")
min_error_location = np.unravel_index(np.argmin(pivot_df.values), pivot_df.shape)
best_delta = pivot_df.index[min_error_location[0]]
best_beta = pivot_df.columns[min_error_location[1]]

square_fig = (7,4)
fig, ax = plt.subplots(1, figsize=square_fig)
contour = ax.contourf(pivot_df.columns, pivot_df.index, pivot_df.values, levels=25, linestyles=None)
#contour.set_clim(2, 25)
fig.colorbar(contour, ax=ax)
label = r'Best score: $\beta$ - ' +f"{best_beta:.2f}, " +r"$\delta$ - " + f"{best_delta:.2f}"
ax.plot(best_beta, best_delta, 'r*', markersize=12, label=label)
ax.set_xlabel(r'TSG/OG parameter - $\beta$')
ax.set_ylabel(r'Acceptance Modulator - $\delta$')
ax.set_title(r'Error Surface for $\beta$ and $\delta$ for Tetraploid Samples')
ax.legend()
fig.savefig(pjoin(base_dir, "img/tet_tsg_and_delta_scan.png"), dpi=300, bbox_inches="tight")
fig.savefig(pjoin(base_dir, "img/tet_tsg_and_delta_scan.pdf"), bbox_inches="tight")


# %%
pivot_df = wgd_0_df.pivot_table(index="delta", columns="beta", values="wd")
min_error_location = np.unravel_index(np.argmin(pivot_df.values), pivot_df.shape)
best_delta = pivot_df.index[min_error_location[0]]
best_beta = pivot_df.columns[min_error_location[1]]

square_fig = (7,4)
fig, ax = plt.subplots(1, figsize=square_fig)
contour = ax.contourf(pivot_df.columns, pivot_df.index, pivot_df.values, levels=25, linestyles=None)
#contour.set_clim(-0.2, 0.2)
fig.colorbar(contour, ax=ax)
label = r'Best score: $\beta$ - ' +f"{best_beta:.2f}, " +r"$\delta$ - " + f"{best_delta:.2f}"
ax.plot(best_beta, best_delta, 'r*', markersize=12, label=label)
ax.set_xlabel(r'TSG/OG parameter - $\beta$')
ax.set_ylabel(r'Acceptance Modulator - $\delta$')
ax.set_title(r'Error Surface for $\beta$ and $\delta$ for Diploid Samples')
ax.legend()
fig.savefig(pjoin(base_dir, "img/dip_tsg_and_delta_scan.png"), dpi=300, bbox_inches="tight")
fig.savefig(pjoin(base_dir, "img/dip_tsg_and_delta_scan.pdf"), bbox_inches="tight")

# %%
# Get the ploidies of observed samples
ploidy_0_df = pd.read_csv(pjoin(base_dir, "data", "wgd_negative_ploidy.tsv"), sep="\t")
ploidy_0_df["ploidy"] = ploidy_0_df["ploidy_major_cn"] + ploidy_0_df["ploidy_minor_cn"]
mean_0 = 1#np.mean(ploidy_0_df.ploidy.values)

ploidy_1_df = pd.read_csv(pjoin(base_dir, "data", "wgd_positive_ploidy.tsv"), sep="\t")
ploidy_1_df["ploidy"] = ploidy_1_df["ploidy_major_cn"] + ploidy_1_df["ploidy_minor_cn"]
mean_1 = 1#np.mean(ploidy_1_df.ploidy.values)

# Take the weighted sum of the two heatmaps
# Ensure the DataFrames are aligned on "beta" and "delta"
merged_df = pd.merge(wgd_0_df[['beta', 'delta', 'wd']], wgd_1_df[['beta', 'delta', 'wd']], on=['beta', 'delta'], suffixes=('_0', '_1'))
merged_df["weighted_wd"] = merged_df["wd_0"]/mean_0 + merged_df["wd_1"]/mean_1

pivot_df = merged_df.pivot_table(index="delta", columns="beta", values="weighted_wd")
min_error_location = np.unravel_index(np.argmin(pivot_df.values), pivot_df.shape)
best_delta = pivot_df.index[min_error_location[0]]
best_beta = pivot_df.columns[min_error_location[1]]

square_fig = (7,4)
fig, ax = plt.subplots(1, figsize=square_fig)
contour = ax.contourf(pivot_df.columns, pivot_df.index, pivot_df.values, levels=25, linestyles=None)
#contour.set_clim(4,12)
fig.colorbar(contour, ax=ax)
label = r'Best score: $\beta$ - ' +f"{best_beta:.2f}, " +r"$\delta$ - " + f"{best_delta:.2f}"
ax.plot(best_beta, best_delta, 'r*', markersize=12, label=label)
ax.set_xlabel(r'TSG/OG parameter - $\beta$')
ax.set_ylabel(r'Acceptance Modulator - $\delta$')
ax.set_title(r'Effect of $\beta$ and $\delta$ for Diploid and Tetraploid Samples')
ax.legend()
fig.savefig(pjoin(base_dir, "img/all_tsg_and_delta_scan.png"), dpi=300, bbox_inches="tight")
fig.savefig(pjoin(base_dir, "img/all_tsg_and_delta_scan.pdf"), bbox_inches="tight")

# CNPs for best weighted_wd
cns_1 = get_cnps(f"{int(best_beta)}_{best_delta}", 1)
cns_1["sample_id"] = "Tet. - " + r"$\beta$: " + f"{int(best_beta)}, " + r"$\delta$: "+ f"{best_delta:.2f}"
cns_0 = get_cnps(f"{int(best_beta)}_{best_delta}", 0)
cns_0["sample_id"] = "Dip. - " + r"$\beta$: " + f"{int(best_beta)}, " + r"$\delta$: "+ f"{best_delta:.2f}"

d = {"Observed WGD+": group_obs_dict[1], "Observed WGD-": group_obs_dict[0], "Tetraploid Synthetic": cns_1, "Diploid Synthetic": cns_0}
combined_df = pd.concat(d.values())
import matplotlib as mpl
colors = ["#66c2a5",  "#8da0cb", "#fc8d62","#e78ac3"]
fig, ax = cns.fig_lines(combined_df, cn_columns="total_cn", colors=colors)
plt.ylim(1,6.5)
ax.legend().remove()
fig.set_size_inches(16,2)
ax.set_title("Copy number profiles on the linear genome, grouped by 3MB bins")
ax.set_xlabel("Chromosomes")
ax.set_ylabel("Total CN")
handles, labels = ax.get_legend_handles_labels()
new_labels = [l for l in d.keys()]
ax.legend(handles, new_labels, ncols=4, loc="upper left")

#fig, ax = cns.fig_lines(combined_df, cn_columns="total_cn", colors=)
plt.savefig(pjoin(base_dir, "img/cnps_tsg_best_joint_dist.png"), dpi=300, bbox_inches="tight")
plt.savefig(pjoin(base_dir, "img/cnps_tsg_best_joint_dist.pdf"), bbox_inches="tight")
