# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
from os.path import join
import json

# %%
n = 1000
simcha_path = "/projects/ag-schwarzr/project-simcha/simcha/SimChA"
template_path = f"{simcha_path}/../configs/wgd_only.json"

def run_simcha(output,params_file):
    cmd = f"dotnet run --no-build --project {simcha_path} -e -R {n} -O {output} -C {params_file} -D {simcha_path}/../data/hg19"
    print(cmd)
    subprocess.run([cmd], shell=True)

p_k0 = 0.7
p_k1 = 0.295
p_k2 = 0.005
real_probs = np.array([p_k0, p_k1, p_k2])
def l2(p, q):
    return np.sqrt(np.sum((p - q)**2))

def update_params(input_file, output_file, stress, mut):
    with open(input_file, 'r') as f:
        data = json.load(f)

    # Update the Stress value
    data["Fitness"]["Stress"] = stress

    # Update the WholeGenomeDoubling Prob value
    data["EvoParams"]["MutationRate"] = mut

    # Write updated data to output file
    with open(output_file, 'w') as f:
        json.dump(data, f, indent=4)

def get_synthetic_counts(events_file, n_samples, n_classes=3):
    counts = np.zeros(n_classes)
    events = pd.read_csv(events_file, sep="\t")
    # Filter events for WholeGenomeDoubling
    wgd_events = events[events["event_type"] == "WholeGenomeDoubling"]
    # Count the number of WGD events per sample_id
    wgd_counts = wgd_events["sample_id"].value_counts()
    zero_count = n_samples - wgd_events.sample_id.nunique()
    # Update counts array
    for count in wgd_counts:
        if count >= len(counts):
            counts[-1] += 1
        else:
            counts[count] += 1
    counts[0] = zero_count
    
    return counts

res = []
i = 0
steps=2
input_file = ""
for mut_rate in np.linspace(0.0001, 0.1, steps):
    for alpha in np.linspace(0.0, 1, steps):
        tmp_file = f"{simcha_path}/../tmp/params_{alpha}_{mut_rate}.json"
        params = update_params(template_path, tmp_file, alpha, mut_rate)
        out = f"{simcha_path}/../scan/{alpha}_{mut_rate}"
        run_simcha(out, params)
        subprocess.run([f"rm {tmp_file}"], shell=True)
        syn_counts = get_synthetic_counts(join(out, "events.tsv"), n)
        syn_freq = syn_counts / sum(syn_counts)
        res.append([mut_rate, alpha, l2(real_probs, syn_freq)])
        i += 1
        print(f"{i/steps**2:.2%}", end='\r')

# %%i
df = pd.DataFrame(res, columns=["mut_rate", "alpha", "error"])
df.to_csv(f"{simcha_path}/../heatmap.png", sep="\t")
pivot_df = df.pivot_table(index="mut_rate", columns="alpha", values="error")

# Find the minimum error value and its corresponding p_wgd and alpha
min_error = pivot_df.min().min()
min_error_location = np.unravel_index(np.argmin(pivot_df.values), pivot_df.shape)
best_p_wgd = pivot_df.index[min_error_location[0]]
best_alpha = pivot_df.columns[min_error_location[1]]

# Plot the error surface
plt.figure(figsize=(10, 8))
contour = plt.contourf(pivot_df.columns, pivot_df.index, pivot_df.values, levels=50, cmap='viridis')
plt.colorbar(contour)
plt.plot(best_alpha, best_p_wgd, 'r*', markersize=12, label='Best score')
plt.xlabel('Alpha')
plt.ylabel('Mutation Rage')
plt.title('Error Surface')
plt.legend()
plt.savefig(f"{simcha_path}/../heatmap.png", dpi=300, bbox_inches="tight")
plt.close()

print(f"Best mutation rate: {best_p_wgd}")
print(f"Best alpha: {best_alpha}")
densities = calc_density(best_p_wgd, best_alpha)
print(f"Density", )
