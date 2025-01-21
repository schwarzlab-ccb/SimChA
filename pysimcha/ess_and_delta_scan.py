#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from os.path import join as pjoin
#import simcha.utils as siu

square_fig = (7,4)

# %%
#p_hemi = 0.184364
#p_null = 0.000308
p_null = 0.0001
#p_rest = 1.0 - p_hemi - p_null
#pcawg_freq = np.array([p_hemi, p_null, p_rest])
#print(pcawg_freq)
def err_func(sim_null_freq, p_null):
    return sim_null_freq - p_null

# %%
#sim_dir = pjoin(siu.get_SimChA_out(), "iscb_results", "results_all_events_noWGD")
"""
sim_dir = "../results_ISMB_ess_and_delta_scan"
subdirs = [d for d in os.listdir(sim_dir) if os.path.isdir(pjoin(sim_dir, d))]
sim_res = []
for i, subdir in enumerate(subdirs):
    print(f"Folders processed: {i}/{len(subdirs)}")
    params = subdir.split("_")
    ess = float(params[0])
    delta = float(params[1])
    if not os.path.exists(pjoin(sim_dir, subdir, "clones.tsv")):
        continue
    clones_df = pd.read_csv(pjoin(sim_dir, subdir, "clones.tsv"), sep="\t")
    hemi_samples = np.mean(clones_df.hemizygosity.values)
    null_samples = np.mean(clones_df.nullizygosity.values)
    rest_samples = 1.0 - hemi_samples - null_samples
    err = err_func(null_samples, p_null)
    #err = err_func(np.array([hemi_samples, null_samples, rest_samples]))
    sim_res.append([ess, delta, err, abs(err), hemi_samples, null_samples, rest_samples, np.mean(clones_df.ploidy.values)])

df = pd.DataFrame(sim_res, columns=["ess", "delta", "error", "abs_error", "hemi", "null", "rest", "ploidy"])
df.to_csv("ess_and_delta_scan.tsv", sep="\t", index=False)
"""

#%%
df = pd.read_csv("ess_and_delta_scan.tsv", sep="\t")
df["error"] = err_func(df["null"], p_null)
print(df.head())
pivot_df = df.pivot_table(index="delta", columns="ess", values="error")
pivot_abs = df.pivot_table(index="delta", columns = "ess", values="abs_error")
min_error_location = np.unravel_index(np.argmin(pivot_abs.values), pivot_abs.shape)
best_delta = pivot_abs.index[min_error_location[0]]
best_alpha = pivot_abs.columns[min_error_location[1]]


# Plot the error surface
fig, ax = plt.subplots(1, figsize=square_fig)
levels = np.linspace(-2*p_null, 2*p_null, 31)
contour = ax.contourf(pivot_df.columns, pivot_df.index, pivot_df.values, levels=levels, cmap="PuOr", linestyles=None, extend="both")
contour.set_clim(levels[0], levels[-1])
cbar = fig.colorbar(contour, ax=ax, extend='both')
cbar.set_label(r"$\mu_{\mathrm{sim}} - \mu_{\mathrm{obs}}$")
label = r'Best score: $\gamma$ - ' +f"{best_alpha:.3f}, " +r"$\delta$ - " + f"{best_delta:.3f}"
ax.plot(best_alpha, best_delta, 'r*', markersize=12, label=label)
import math as m
x = np.linspace(0, 1000, 10000)
y = [(1500-i)**0.38/2.5 for i in x]
#y = [m.log(900-i,2)/2 for i in x]
#y = [m.tan((450-i)/450 * m.pi/2)*2.6 for i in x]
#ax.plot(x, y, color='k', ls="--", alpha=0.6, label=f"y=-x/200+5")
ax.set_xlabel(r'Essentiality - $\gamma$')
ax.set_ylabel(r'$\delta$')
ax.set_ylim(0,4.8)
ax.legend(loc="lower left")
ax.set_title(r"Essentiality and $\delta$-parameter's effect on WGD- Nullizygosity")
fig.savefig("../img/nonWGD_ess_and_delta_scan.png", dpi=300, bbox_inches="tight")
fig.savefig("../img/nonWGD_ess_and_delta_scan.pdf")
#siu.save_fig_out("nonWGD_ess_scan", fig)
plt.show()
#%%
#print(f"PCAWG:\t hemi: {pcawg_freq[0]:6f}; nulli: {pcawg_freq[1]:6f}; rest: {pcawg_freq[2]:6f}")
"""
print(f"PCAWG:\t nulli: {p_null:6f} ")
best_error = df.query("error == error.min()")
print(f"Best:\t nulli: {best_error.null.values[0]:.6f}; hemi: {best_error.hemi.values[0]:.6f}; rest: {best_error.rest.values[0]:.6f}")
worst_error = df.query("error == error.max()")
print(f"Worst:\t nulli: {worst_error.null.values[0]:.6f}; hemi: {worst_error.hemi.values[0]:.6f}; rest: {worst_error.rest.values[0]:.6f}")
"""
# %%
"""
fig, ax = plt.subplots(1, figsize=square_fig)
ax.scatter(df["ess"], df["ploidy"], alpha=1)
ax.set_xlabel(r"Essentiality parameter - $\gamma$")
ax.set_ylabel("Ploidy")
ax.set_title(r"Effect of Essentiality on WGD- ploidies")
fig.savefig("../img/nonWGD_ess_and_delta_scan_ploidy.png", dpi=300, bbox_inches="tight")
fig.savefig("../img/nonWGD_ess_and_delta_scan_ploidy.pdf")
#siu.save_fig_out("nonWGD_ess_scan_ploidy", fig)
plt.show()
"""
