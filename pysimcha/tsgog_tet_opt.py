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

status = 1

group_obs_dict[status] = pd.read_csv(pjoin(base_dir, "data",f"grouped_obs_cnps_wgd_{status}.tsv"), sep="\t")
mean_total_cn = group_obs_dict[status].groupby('chrom')['total_cn'].transform('mean')
group_obs_dict[status]["diff_mean"] = group_obs_dict[status]['total_cn'] - mean_total_cn

def wasserstein_dist(obs_df, sim_df):
    wasserstein_dist = 0
    for chrom in obs_df.chrom.unique():
        chrom_df = sim_df.query(f"chrom == '{chrom}'")
        wasserstein_dist += wasserstein_distance(
            obs_df.query(f"chrom == '{chrom}'").diff_mean,
            sim_df.query(f"chrom == '{chrom}'").diff_mean) if len(chrom_df) > 0 else float("inf")
    return wasserstein_dist

out_name = pjoin(base_dir, "out", "tsgog_tet_results_mean.tsv")
if os.path.exists(out_name):
    results_df = pd.read_csv(out_name, sep="\t")
else:
    results_df = pd.DataFrame(columns=["beta", "delta", "status", "wd"])
cns_dict = {}
groups_df = {}
for _, curr_dir in enumerate([tetraploid_base_dir]):
    groups_df[status] = {}
    cns_dict[status] = {}
    for i, subdir in enumerate(subdir_dict[status]):
        print(i)
        params = subdir.split("_")
        beta = float(subdir.split("_")[0])
        delta = float(subdir.split("_")[1])
        if len(results_df.query("beta == @params[0] and delta == @params[1]")) > 0:
            continue
        current_dir = pjoin(curr_dir, subdir)
        if not os.path.exists(pjoin(current_dir, "cns_3MB.tsv")):
            continue
        cns_dict[status][subdir] = cns.load_cns(pjoin(current_dir, "cns_3MB.tsv"))
        cns_dict[status][subdir].query("chrom != 'chrX' and chrom != 'chrY'", inplace=True)
        groups_df[status][subdir] = cns.group_samples(cns_dict[status][subdir])
        groups_df[status][subdir]["total_cn"] = (
            groups_df[status][subdir]["hap1_cn"] 
            + groups_df[status][subdir]["hap2_cn"])
        groups_df[status][subdir]["sample_id"] = "Synthetic Tetraploid" if status == 1 else "Synthetic Diploid"
        mean_total_cn = groups_df[status][subdir].groupby('chrom')['total_cn'].transform('mean')
        groups_df[status][subdir]["diff_mean"] = groups_df[status][subdir]["total_cn"] - mean_total_cn
        dist = wasserstein_dist(group_obs_dict[status], groups_df[status][subdir])
        print(dist)
        results_df = pd.concat([results_df, pd.DataFrame({"beta": [beta], "delta": [delta], "status": [status], "wd": [dist]})])
results_df.to_csv(out_name, sep="\t", index=False)
