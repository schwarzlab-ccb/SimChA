# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import cns.data_utils as cdu
#import cns
import os
from os.path import join
import subprocess

#%%
"""samples_df, cns_df = cdu.main_load("imp", use_filter = False)
samples_df.head()
# %%

ploidy_df = cns.main_ploidy(cns_df)
ploidy_df["ploidy"] = ploidy_df["ploidy_major_cn"] + ploidy_df["ploidy_minor_cn"]
ploidy_df.to_csv("ploidy_all_samples.tsv", sep="\t", index=False)
"""
# %%
ploidy_df = pd.read_csv("ploidy_all_samples.tsv", sep="\t")
bins = np.linspace(0, 10, 100)
obs_dist, bin_edges = np.histogram(ploidy_df.ploidy.values, bins=bins, density=False)
# Add a pseudo-count to each bin
obs_dist += 1
obs_probs = obs_dist/obs_dist.sum()
# %%
def err_func(mean_ploidy):
    return (mean_ploidy - 2)
def log_likelihood(simulated, observed_probs, bin_edges):
    sim_count, _ = np.histogram(simulated, bins=bin_edges, density=False)
    sim_count += 1
    sim_probs = sim_count/sim_count.sum()
    ll = np.sum(observed_probs * np.log(sim_probs))
    return -ll

#%%
obs_ratios = pd.Series({"Arm": 0.980, "CentromereBound": 1.518, "Internal": 2.217, "Tail": 0.946, "Chrom": 0.420})
events = ["Arm", "CentromereBound", "Internal", "Tail", "Chrom"]
# I want to compare the results of the two methods
def ratio_dist(event_ratios, obs_ratios):
    return np.sum(np.square(event_ratios - obs_ratios))

def get_event_ratios(events_df):
    event_ratio = {}
    total_losses = 0
    total_gains = 0
    wgd_samples = events_df.query("event_type == 'WholeGenomeDoubling'").sample_id.unique()
    events_df = events_df[events_df.sample_id.isin(wgd_samples)]
    """
    event_counts = events_df['event_type'].value_counts()
    for event in events:
        losses = event_counts.get(f"{event}Deletion", 0)
        gains = event_counts.get(f"{event}Duplication", 0)
    """
    for event in events:
        losses = len(events_df.query(f"event_type == '{event}Deletion'"))
        gains = len(events_df.query(f"event_type == '{event}Duplication'"))
        event_ratio[event] = gains / losses if losses != 0 else float('inf')
        total_losses += losses
        total_gains += gains
    
    return event_ratio

# %%
"""
sim_dir = "../results_event_cost"
subdirs = [d for d in os.listdir(sim_dir) if os.path.isdir(join(sim_dir, d))]
res = []
import json
for i, subdir in enumerate(subdirs):
    print(f"Folders processed: {i}/{len(subdirs)}")
    params = subdir.split("_")
    event_cost = float(params[0])
    alpha = float(params[1])
    clones_file = join(sim_dir, subdir, "clones.tsv")
    sim_params = join(sim_dir, subdir, "sim_params.json")
    with open(sim_params, "r") as f:
        file = json.load(f)
    p_wgd = file['Signatures']['CNVs']['Events'][-1]['Prob']
    if not os.path.exists(clones_file):
        continue
    ploidy = pd.read_csv(clones_file, sep="\t").ploidy.values
    ll = log_likelihood(ploidy, obs_probs, bin_edges)
    events_df = pd.read_csv(join(sim_dir, subdir, "events.tsv"), sep="\t")
    event_ratios = get_event_ratios(events_df)
    d_ratios = ratio_dist(pd.Series(event_ratios), obs_ratios)
    print(d_ratios)
    res.append([event_cost, alpha, ll, p_wgd, d_ratios])
"""
# %%
square_fig = (7,4)
#df = pd.DataFrame(res, columns=["event_cost", "alpha", "ll", "p_wgd", "ratios"])
#df["total"] = df["ll"] + df["ratios"] 
#df.to_csv("event_cost_0_10_with_ratios.tsv", sep="\t", index=False)
df = pd.read_csv("event_cost_0_10_with_ratios.tsv",sep="\t")
df["log"] = df["total"]
pivot_df = df.pivot_table(index="event_cost", columns="alpha", values="log")
# Find the minimum error value and its corresponding p_wgd and alpha
min_error_location = np.unravel_index(np.argmin(pivot_df.values), pivot_df.shape)
best_r_loss = pivot_df.index[min_error_location[0]]
best_alpha = pivot_df.columns[min_error_location[1]]
best_pwgd = df[(df['alpha'] == best_alpha) & (df['event_cost'] == best_r_loss)]['p_wgd'].values[0]

# Find the rows with the min z values
min_z_rows = df.loc[df.groupby("alpha")["total"].idxmin()]
x = min_z_rows["alpha"].values.reshape(-1, 1)
y = min_z_rows["event_cost"].values

from sklearn.linear_model import LinearRegression
model = LinearRegression()
model.fit(x, y)
m = model.coef_[0]
c = model.intercept_

# Plot the error surface
fig, ax = plt.subplots(1, figsize=square_fig)
levels = np.linspace(3.6, 10, 31)
contour = ax.contourf(pivot_df.columns, pivot_df.index, pivot_df.values, levels=levels, linestyles=None, extend='max')
contour.set_clim(levels[0], levels[-1])
cbar = fig.colorbar(contour, ax=ax, extend='max')
label = r'Best score: $\alpha$ - ' +f"{best_alpha:.3f}, " + "m - " + f"{best_r_loss:.3f}, " r"$p_{\mathrm{wgd}}$ - " + f"{best_pwgd:.3f}"
ax.plot(best_alpha, best_r_loss, 'r*', markersize=12, label=label)

x_vals = np.linspace(df["alpha"].min(), df["alpha"].max(), 500)
y_vals = m * x_vals + c
ax.plot(x_vals, y_vals, color="w", label=f"y = {m:.3f}x+ {c:.3f}")

ax.set_xlabel(r'Stress parameter - $\alpha$')
ax.set_ylabel(r'Multiplication Factor - $m$')
ax.set_title(r'Effect of $\alpha$ and $m$ on WGD+ Samples')
ax.legend()
fig.savefig("../img/event_cost_total.png", dpi=300, bbox_inches="tight")
fig.savefig("../img/event_cost_total.pdf")
"""
pivot_df = df.pivot_table(index="event_cost", columns="alpha", values="ratios")
# Find the minimum error value and its corresponding p_wgd and alpha
min_error_location = np.unravel_index(np.argmin(pivot_df.values), pivot_df.shape)
best_r_loss = pivot_df.index[min_error_location[0]]
best_alpha = pivot_df.columns[min_error_location[1]]

# Plot the error surface
fig, ax = plt.subplots(1, figsize=square_fig)
contour = ax.contourf(pivot_df.columns, pivot_df.index, pivot_df.values, levels=100, cmap='PuOr', linestyles=None)
fig.colorbar(contour, ax=ax)
label = r'Best score: $\alpha$ - ' +f"{best_alpha:.3f}, " +r"$R_{\mathrm{loss}}$ - " + f"{best_r_loss:.3f}"
ax.plot(best_alpha, best_r_loss, 'r*', markersize=12, label=label)
ax.set_xlabel(r'Stress parameter - $\alpha$')
ax.set_ylabel(r'Event Cost')
ax.set_title(r'Stress vs Post-WGD Event Cost Error Surface')
ax.legend()
fig.savefig("../img/event_cost_ratios_long.png", dpi=300, bbox_inches="tight")
fig.savefig("../img/event_cost_ratios_long.pdf")
"""
