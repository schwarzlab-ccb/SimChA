import numpy as np
from numba import njit
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution
import json
from os.path import join
import subprocess
from multiprocessing import cpu_count
import uuid
import shutil
import os

import psutil
print(f"SLURM CPU cores allocated: {psutil.cpu_count(logical=False)}")

num_cores = 90
output_file = "joint_optimizer_results_no_tails_haplo_suf.tsv"


data_dir = "../out/pcawg"
pcawg_clones = pd.read_csv(join(data_dir, "clones.tsv"), sep="\t")
ploidy_min = 1.2
ploidy_max = 9

ploidy_bins = np.linspace(ploidy_min, ploidy_max, 50)
ploidy_values = pcawg_clones.ploidy.values

wgd_bins = np.linspace(0, 4, 4)
wgd_values = np.zeros(len(ploidy_values))

H, x_edges, y_edges = np.histogram2d(ploidy_values, wgd_values, bins=[ploidy_bins, wgd_bins])
H += 1
joint_obs = H/np.sum(H)
joint_obs_df = pd.DataFrame(joint_obs, index=x_edges[:-1], columns = y_edges[:-1]) 
joint_obs_df.to_csv("joint_pcawg_df.tsv", sep="\t", index=False)

def log_likelihood(joint_sim, joint_obs):
    # Flatten then arrays
    obs_flat = joint_obs.flatten()
    sim_flat = joint_sim.flatten()

    sim_flat = np.maximum(sim_flat, 1e-10)

    ll = np.sum(obs_flat * np.log(sim_flat))
    return -ll

simcha_path = "/projects/ag-schwarzr/project-simcha/simcha/"
config_file = join(simcha_path, "configs", "wgd_and_loss.json")

n_events = 9
theta = 1
c_wgd = 0.35
def get_wgd_weight(alpha, delta):
    return np.exp(delta*alpha)*c_wgd/n_events

def run_simcha(params):
    
    print(f"Worker PID: {os.getpid()}, Params: {params}")
    alpha, r_loss, delta = params
    t = 2*r_loss + 2
    r_wgd = get_wgd_weight(alpha, delta)
    wgd = r_wgd*t/(1-r_wgd)#get_wgd_weight(alpha, theta)
    print(f"new run: {alpha}, {r_loss}, {delta}, {wgd}")
    with open(config_file, 'r') as f:
        config = json.load(f)
    config['Fitness']['Stress'] = alpha
    config['Fitness']['TsgOg'] = 0
    config['Fitness']['Essentiality'] = 1.0-alpha
    config['Fitness']['TotalStrength'] = delta
    config['Fitness']['Haploinsufficiency'] = False
    config['Signatures']['CNVs']['Events'][0]['Prob'] = r_loss
    config['Signatures']['CNVs']['Events'][1]['Prob'] = r_loss
    #config['Signatures']['CNVs']['Events'][2]['Prob'] = r_loss
    config['Signatures']['CNVs']['Events'][-1]['Prob'] = wgd
    config['EvoParams']['WithFitness'] = True
    config['EvoParams']['EvolveInTime'] = False
    config['EvoParams']['ThetaFitness'] = theta
    config['EvoParams']['MutationRate'] = 1
    pid = str(uuid.uuid4())
    run_config = join(simcha_path, "scripts", "temp", f"config_{pid}.json")
    with open(run_config, 'w') as f:
        json.dump(config, f, indent = 4)
    cmd = f"dotnet run --no-build --project {simcha_path}/SimChA -C {run_config} -D {simcha_path}/data/hg19 -e -R 2000 -O {join(simcha_path, 'scripts', 'temp', pid)} --light"
    #proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #proc.wait()
    subprocess.run([cmd], shell=True, check=True)
    return pid
    
max_wgd_weight = 0.25
def alpha_bound(delta):
    # limit alpha such that the weight of wgd isn't larger than max_wgd_weight
    return np.log(max_wgd_weight*n_events/c_wgd)/delta

def get_wgd_status(events_df):
    n_samples = events_df.sample_id.nunique()
    # Get the samples with WGD
    wgd_events = events_df.query("event_type == 'WholeGenomeDoubling'")
    nonwgd_counts = np.zeros(n_samples - wgd_events.sample_id.nunique())
    wgd_counts = np.array(wgd_events.sample_id.value_counts().values)
    all_counts = np.append(nonwgd_counts, wgd_counts)
    return all_counts

def err_func(params):
    alpha, r_loss, delta = params
    # Dynamic adjustment of bounds for alpha based on the delta param
    max_alpha = alpha_bound(delta)
    if not alpha <= max_alpha:
        return float('inf')
    pid = run_simcha(params)
    try:
        ploidy_sim = pd.read_csv(join("temp", pid, "clones.tsv"), sep="\t").ploidy.values
        events_df = pd.read_csv(join("temp", pid, "events.tsv"), sep="\t")
        wgd_sim = get_wgd_status(events_df)
        H, x_edges, y_edges = np.histogram2d(ploidy_sim, wgd_sim, bins=[ploidy_bins, wgd_bins])
        H += 1
        joint_sim = H/np.sum(H)
        return log_likelihood(joint_obs, joint_sim)
    finally:
        # Safe cleanup of temp files
        shutil.rmtree(join("temp", pid), ignore_errors=True)
        os.remove(join("temp", f"config_{pid}.json"))

initial_guess = [0.05, 2.5, 1]

bounds = [(0, alpha_bound(1)), (0, 10), (0.25, 40)]

if __name__ == "__main__":
    result = differential_evolution(err_func, bounds, x0=initial_guess, workers=num_cores)

    if result.success:
        print("Optimizer succeeded!")
        best_alpha, best_r_loss, best_delta = result.x
        dict_results = [[best_alpha, best_r_loss, best_delta]]
        df = pd.DataFrame(dict_results, columns=["alpha", "r_loss", "delta"])
        df.to_csv(output_file, sep="\t")
    else:
        print("Optimizer failed: ", result.message)
