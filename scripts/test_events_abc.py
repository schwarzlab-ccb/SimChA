#!/usr/bin/env python3
import subprocess
import numpy as np
import os
import argparse
import platform
from datetime import datetime
import uuid
import json
import matplotlib.pyplot as plt
from pyabc import ABCSMC, RV, Distribution, settings, visualization, sampler
from pyabc.visualization import plot_kde_matrix
from pyabc.populationstrategy import ConstantPopulationSize
import pyabc

# This python program is designed to test if the pyABC sampling process can infer
# the parameter set used by SimChA to generate a simulated cohort.

settings.set_figure_params('pyabc')

def update_params_file(params):
    # Create the temporary parameter file with a random parameter selected as the suffix
    foldername = f"{uuid.uuid4()}" 
    path = f"{pwd}/temp/{foldername}"
    subprocess.run([f"mkdir -p {path}"], shell = True)
    subprocess.run([f"cp simple_params.json {path}"], shell = True)
    # Create the output subdirectory
    subprocess.run([f"mkdir -p {path}/out"], shell=True)
    file_path = f"{path}/simple_params.json"

    # Update the parameter file
    with open(file_path, 'r', encoding='utf-8') as json_file:
        configs = json.load(json_file)

    # Events are unfortunately an array in the parameters file
    # Also round the weights to the nearest 3 decimal places
    ndp = 3
    # The length-scales have to be rounded to integer status 
    # ChromDeletion
    configs["Signatures"]["CNs"]["Events"][0]["Prob"] = round(float(params["w_chrom_del"]), ndp)
    # ChromDuplication
    configs["Signatures"]["CNs"]["Events"][1]["Prob"] = round(float(params["w_chrom_dup"]), ndp)
    # InternalDuplication
    configs["Signatures"]["CNs"]["Events"][2]["Prob"] = round(float(params["w_int_dup"]), ndp)
    configs["Signatures"]["CNs"]["Events"][2]["Size"] = int(params["l_int_dup"]*100_000)
    # InternalDeletion
    configs["Signatures"]["CNs"]["Events"][3]["Prob"] = round(float(params["w_int_del"]), ndp)
    configs["Signatures"]["CNs"]["Events"][3]["Size"] = int(params["l_int_del"]*100_000)
    # InternalInversion
    configs["Signatures"]["CNs"]["Events"][4]["Prob"] = round(float(params["w_int_inv"]), ndp)
    configs["Signatures"]["CNs"]["Events"][4]["Size"] = int(params["l_int_inv"]*100_000)
    # InvertedDuplication
    configs["Signatures"]["CNs"]["Events"][5]["Prob"] = round(float(params["w_inv_dup"]), ndp)
    configs["Signatures"]["CNs"]["Events"][5]["Size"] = int(params["l_inv_dup"]*100_000)
    # BreakageFusionBridge    
    configs["Signatures"]["CNs"]["Events"][6]["Prob"] = round(float(params["w_bfb"]), ndp)
    # TailDeletion
    configs["Signatures"]["CNs"]["Events"][7]["Prob"] = round(float(params["w_tail_del"]), ndp)
    # Translocation
    configs["Signatures"]["CNs"]["Events"][8]["Prob"] = round(float(params["w_transloc"]), ndp)
    configs["Signatures"]["CNs"]["Events"][8]["Size"] = int(params["l_transloc"]*100_000)
    # Whole-Genome Doubling
    configs["Signatures"]["CNs"]["Events"][9]["Prob"]  = float(params["w_wgd"])
    
    with open(file_path, 'w', encoding="utf-8") as json_file:
        json.dump(configs, json_file)
    # Return the path to the config file
    return path


def run_simcha(params, genes_path, cohort_path, repeats):
    param_file_path = update_params_file(params)

    cmd = f"dotnet run --no-build --project SimChA -C {param_file_path}/simple_params.json -R {repeats} -O {param_file_path}/out --optimization events -D {genes_path} -P {cohort_path} --autosomes-only"
    output = subprocess.check_output([cmd], universal_newlines=True, shell=True)
    # SimChA produces as its output the Euclidean sum of Wasserstein distances for each of the 
    # characteristic features of cancer genomes, printing the double to the command 
    last_line = output.strip().split("\n")[-1]

    # Delete the temporary folder and files
    subprocess.run([f"rm -rf {param_file_path}"], shell=True)
    # Return the distance SimChA calculated
    return float(last_line.split(":")[1].strip())

# Distance function required by pyABC, but note that SimChA does the actual calculation, and the observed data is set to distance 0.
def distance(x,y):
    return abs(x["distance"] - y["distance"])

# Function to generate the ground-truth SimChA simulated dataset
def generate_cohort(param_file, genes_path, out_path, repeats):
    cmd = f"dotnet run --no-build --project SimChA -C {param_file} -R {repeats} -O {out_path} -D {genes_path}"
    subprocess.run([cmd], shell=True)
    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="pyABC program to fit parameters in SimChA ")
    #parser.add_argument('-R', "--repeats", type=int, default=500, help="Number of SimChA simulated samples to generate for each pyABC sample")
    parser.add_argument('-N', "--name",type=str, default="test_events_abc_results", help="Name for output directory to put SQL database produced by pyABC and the posterior plot produced")
    parser.add_argument("-r", "--repeats", type=int, default=2000, help="Number of samples/repeats to generate for each SimChA run")
    parser.add_argument('--n_procs', type=int, default = 8, help="Number of parallel threads that pyABC will use")
    parser.add_argument("--n_pop", type=int, default=150, help="Population size, i.e. number of accepted samples (particles) to move on to new generation")
    parser.add_argumetn("--max_gen", type=int, default=10, help="Maximum number of generations to consider")
    parser.add_argument("--min_eps", type=float, default=0.01, help="Minimum acceptance threshold before ending the ABC Sampling prematurely")
    args = parser.parse_args()

    genes_path = "data/hg19_1000"

    # Build the program once at the beginning in case any changes have been made to SimChA
    subprocess.run(["dotnet build"], shell=True)
 
    # Wrapper function for model so that we can run SimChA with the input dataset and any modified hyperparameters (like number of SimChA samples)
    def model_wrapper(params):
        return {"distance": run_simcha(params, genes_path, cohort_path, args.repeats)}
    
    # Create the output directory for the final plots as well as the SQL database
    pwd = os.getcwd()
    out_dir = args.name
    subprocess.run([f"mkdir -p {out_dir}"], shell=True)

    param_file = f"{pwd}/simple_params.json"
    cohort_dir_path = "out/ground_truth"
    # Generate the simulated data as the ground truth
    generate_cohort(param_file, genes_path, cohort_dir_path, repeats)
    cohort_path = "out/ground_truth/copynumbers.tsv"

    # Uniform prior distributions for the various different properties of the simple events
    limits = dict(
            w_chrom_del = (1, 11),
            w_chrom_dup = (1, 11),
            w_int_dup   = (1, 101),
            w_int_del   = (1, 101),
            w_int_inv   = (1, 101),
            w_inv_dup   = (1, 101),
            w_bfb       = (1, 11),
            w_tail_del  = (1, 11),
            w_wgd       = (0.01, 0.11),
            w_transloc  = (1, 26),
            l_int_dup   = (1, 51), 
            l_inv_dup   = (1, 51),
            l_int_del   = (1, 51),
            l_int_inv   = (1, 51),
            l_transloc  = (1, 101))
    # The length scale of events in is 100kb.
    # Priors are simply a uniform distribution between their two limits
    prior = Distribution(**{key: RV("uniform", a, b-a) for key, (a, b) in limits.items()})
    
    # SimChA calculates the Euclidean-summed Wasserstein distance, so we don't need an observed distance
    observed_data = {"distance": 0.0}
    sampler = sampler.MulticoreEvalParallelSampler(n_procs=args.n_procs)
    abc = ABCSMC(model_wrapper, prior, distance_function=distance, population_size = args.n_pop, sampler = sampler)
    # ABC-SMC output is a SQL database
    db_path = f"{out_dir}/parameters.db"
    abc.new("sqlite:///"+db_path, observed_data)
    
    # Main part of the ABC Sampling
    history = abc.run(minimum_epsilon = args.min_eps, max_nr_populations = args.max_gen)
    
    # ID 
    fig, ax = plt.subplots()
    for t in range(history.max_t + 1):
        df, w = history.get_distribution(t=t)
        visualization.plot_kde_1d(
            df,
            w,
            xmin=0.01,
            xmax=0.11,
            x="w_wgd",
            xname=r"Event Count",
            ax=ax,
            label=f"PDF t={t}",
        )
    ax.legend()
    plt.savefig(f"{out_dir}/posterior_event_count.png")
    df, w = history.get_distribution(m=0, t=history.max_t)
    # the KDE matrix is an array of axes
    arr_ax = plot_kde_matrix(df, w, limits=limits)
    fig = arr_ax[0, 0].get_figure()
    fig.savefig(f"{out_dir}/kde_matrix.png", dpi=300, bbox_inches="tight")

