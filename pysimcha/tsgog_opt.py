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

tetraploid_base_dir = pjoin("..", "results_ISMB_tsg_tetraploid_and_delta_scan")
diploid_base_dir = pjoin("..", "results_ISMB_tsg_diploid_and_delta_scan")

subdir_dict = {}
for status, curr_dir in enumerate([diploid_base_dir, tetraploid_base_dir]):
    subdirs = [d for d in os.listdir(curr_dir) if os.path.isdir(pjoin(curr_dir, d))]
    subdir_dict[status] = subdirs

group_obs_dict = {}

for status in wgd_status:
    group_obs_dict[status] = pd.read_csv(pjoin("..", "data",f"grouped_obs_cnps_wgd_{status}.tsv"), sep="\t")

def wasserstein_dist(obs_df, sim_df):
    wasserstein_dist = 0
    for chrom in obs_df.chrom.unique():
        chrom_df = sim_df.query(f"chrom == '{chrom}'")
        wasserstein_dist += wasserstein_distance(
            obs_df.query(f"chrom == '{chrom}'").total_cn,
            sim_df.query(f"chrom == '{chrom}'").total_cn) if len(chrom_df) > 0 else float("inf")
    return wasserstein_dist

results_df = pd.DataFrame(columns=["beta", "delta", "status", "wd"])
cns_dict = {}
groups_df = {}
for status, curr_dir in enumerate([diploid_base_dir, tetraploid_base_dir]):
    groups_df[status] = {}
    cns_dict[status] = {}
    for i, subdir in enumerate(subdir_dict[status]):
        print(i)
        beta = float(subdir.split("_")[0])
        delta = float(subdir.split("_")[1])
        current_dir = pjoin(curr_dir, subdir)
        cns_dict[status][subdir] = cns.load_cns(pjoin(current_dir, "cns_3MB.tsv"))
        cns_dict[status][subdir].query("chrom != 'chrX' and chrom != 'chrY'", inplace=True)
        groups_df[status][subdir] = cns.group_samples(cns_dict[status][subdir])
        groups_df[status][subdir]["total_cn"] = (
            groups_df[status][subdir]["hap1_cn"] 
            + groups_df[status][subdir]["hap2_cn"])
        groups_df[status][subdir]["sample_id"] = "Synthetic Tetraploid" if status == 1 else "Synthetic Diploid"
        dist = wasserstein_dist(group_obs_dict[status], groups_df[status][subdir])
        print(dist)
        results_df = pd.concat([results_df, pd.DataFrame({"beta": [beta], "delta": [delta], "status": [status], "wd": [dist]})])
results_df.to_csv("dip_tsg_og_optimization_results.tsv", sep="\t", index=False)

# %%
def get_cnps(subdir, status):
    curr_dir = tetraploid_base_dir if status == 1 else diploid_base_dir
    df = cns.load_cns(pjoin(curr_dir, subdir, "cns_3MB.tsv"))
    df.query("chrom != 'chrX' and chrom != 'chrY'", inplace=True)
    cnps = cns.group_samples(df)
    cnps["total_cn"] = cnps["hap1_cn"] + cnps["hap2_cn"]
    return cnps

# %%
df_1 = pd.read_csv("tet_tsg_og_optimization_results.tsv", sep="\t")
df_0 = pd.read_csv("dip_tsg_og_optimization_results.tsv", sep="\t")
results_df = pd.concat([df_1, df_0])#pd.read_csv("tsg_og_optimization_results.tsv", sep="\t")
wgd_1_df = results_df.query("status == 1")
wgd_0_df = results_df.query("status == 0")

best_wgd_1_result = wgd_1_df.query("wd == wd.min()")
beta_1  = best_wgd_1_result.beta.values[0]
delta_1 = best_wgd_1_result.delta.values[0]
#cns_1 = get_cnps(f"{beta_1}_{delta_1}", 1)
cns_1 = groups_df[1][f"{beta_1}_{delta_1}"]
cns_1["sample_id"] = "Tet. - " + r"$\beta$: " + f"{beta_1:.2f}, " + r"$\delta$: "+ f"{delta_1:.2f}"
best_wgd_0_result = wgd_0_df.query("wd == wd.min()")
beta_0  = best_wgd_0_result.beta.values[0]
delta_0 = best_wgd_0_result.delta.values[0]
#cns_0 = get_cnps(f"{beta_0}_{delta_0}", 0)
cns_0 = groups_df[0][f"{beta_0}_{delta_0}"]
cns_0["sample_id"] = "Dip. - " + r"$\beta$: " + f"{beta_0:.2f}, " + r"$\delta$: "+ f"{delta_0:.2f}"

d = {"Obs WGD+": group_obs_dict[1], "Obs WGD-": group_obs_dict[0], "Tet. Syn.": cns_1, "Diploid Syn.": cns_0}
combined_df = pd.concat(d.values())
cns.fig_lines(combined_df, cn_columns="total_cn")
plt.savefig("../img/cnps_tsg_best.png", dpi=300, bbox_inches="tight")
plt.savefig("../img/cnps_tsg_best.pdf", bbox_inches="tight")

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
ax.set_ylabel(r'$\delta$')
ax.set_title(r'Error Surface for $\beta$ and $\delta$ for Tetraploid Samples')
ax.legend()
fig.savefig("../img/tet_tsg_and_delta_scan.png", dpi=300, bbox_inches="tight")
fig.savefig("../img/tet_tsg_and_delta_scan.pdf", bbox_inches="tight")


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
ax.set_ylabel(r'$\delta$')
ax.set_title(r'Error Surface for $\beta$ and $\delta$ for Diploid Samples')
ax.legend()
fig.savefig("../img/dip_tsg_and_delta_scan.png", dpi=300, bbox_inches="tight")
fig.savefig("../img/dip_tsg_and_delta_scan.pdf", bbox_inches="tight")



