#!/usr/bin/env python3
import subprocess
import numpy as np
import os
import argparse
import platform
import datetime as dt
import uuid
import json
import matplotlib.pyplot as plt
from pyabc import ABCSMC, RV, Distribution, settings, visualization, sampler
from pyabc.visualization import plot_kde_matrix
from pyabc.populationstrategy import ConstantPopulationSize
import pyabc

######################
#
# This python program perfoms the pyABC sampling for the fitness parameters used in SimChA
#
######################

######
# Files needed to infer the fitness parameters:
# genes list (here: data/hg19_1000)
# params file (here: ./fitted_params.json)
# Samples already run through SimChA to calculate their fitness components (here: out/PCAWG/clones.tsv)
# Copy-number profiles (here: pcawg_filtered_95_pc.tsv)

settings.set_figure_params('pyabc')

def update_params_file(params):
    # Create the temporary parameter folder and file
    foldername = f"{uuid.uuid4()}"
    path = f"{pwd}/temp/{foldername}"
    subprocess.run([f"mkdir -p {path}"], shell = True)
    subprocess.run([f"cp fitted_params.json {path}"], shell = True)
    # Create the output subdirectory
    subprocess.run([f"mkdir -p {path}/out"], shell=True)
    file_path = f"{path}/fitted_params.json"

    # Update the parameter file
    with open(file_path, 'r', encoding='utf-8') as json_file:
        configs = json.load(json_file)

    # The stress, TSG/OG, and essentiality weights sum to 1, so we have to give the essentiality
    # score a default then normalize the three
    gamma = 0.1
    fitness_sum = float(params["alpha"]) + float(params["beta"]) + gamma

    configs["Fitness"]["Stress"] = float(params["alpha"]) / fitness_sum
    configs["Fitness"]["TsgOg"]  = float(params["beta"]) / fitness_sum
    configs["Fitness"]["Essentiality"] = gamma / fitness_sum
    configs["Fitness"]["TotalStrength"] = float(params["w_strength"])

    with open(file_path, 'w', encoding="utf-8") as json_file:
        json.dump(configs, json_file)
    # Return the path to the config file
    return path

def run_simcha(params, genes_path, cohort_path, bootstrap_path, binned_path, repeats, all_chromosomes):
    param_file_path = update_params_file(params)

    cmd = f"dotnet run --no-build --project SimChA -C {param_file_path}/fitted_params.json -R {repeats} -O {param_file_path}/out --optimization fitness -D {genes_path} -M -B {bootstrap_path} -P {cohort_path} --binned-samples {binned_path}"
    if not all_chromosomes:
        cmd += " --autosomes_only"

    # SimChA produces as its output the Euclidean sum of Wasserstein distances for each of the 
    # characteristic features of cancer genomes, printing the double to the command 
    output = subprocess.check_output([cmd], universal_newlines = True, shell = True)
    last_line = output.strip().split("\n")[-1]

    # Delete the temporary folder and files
    subprocess.run([f"rm -rf {param_file_path}"], shell = True)
    # Return the distance SimChA calculated
    return float(last_line.split(":")[1].strip())

def distance(x,y):
    return abs(x["distance"] - y["distance"])
    
def check_inputs(in_path, genes_path, cohort_path):

    clones_file = "clones.tsv"
    binned_file = "binned_CNs.tsv"

    if !os.path.isdir(in_path):
        subprocess.run([f"mkdir -p in_path"], shell = True)
    # bootstrap file needed for sampling the fitnesses of the mcmc produces clones
    if !os.path.isfile(os.path.join(in_path, clones_file)):
        cmd = f"dotnet run --no-build --project SimChA -P {cohort_path} -O {in_path} -D {genes_path}"
        subprocess.run([cmd], shell = True)
    # To avoid having to generate the binned copy-number profiles for each run of SimChA, we do it once here if the file does not exist.
    if !os.path.isfile(os.path.join(in_path, binned_file)):
        cmd = f"dotnet run --no-build --project SimChA -P {cohort_path} -O {in_path} --bin-samples"
        subprocess.run([cmd], shell = True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="pyABC program to fit parameters in SimChA ")
    #parser.add_argument('-R', "--repeats", type=int, default=500, help="Number of SimChA simulated samples to generate for each pyABC sample")
    parser.add_argument('-n', "--name",type=str, default="fitness_abc_results", help="Name for output directory to put SQL database produced by pyABC and the posterior plot produced")
    parser.add_argument("-w", "--weight", type=float, default=1, help="Total weight associated with the initial guess for the fitness parameters. Higher weight means lower variance, i.e. sample closer to initial guesses.")
    parser.add_argument("-r", "--repeats", type=int, default=250, help="Number of samples to generate for each SimChA run")
    parser.add_argument('--n_procs', type=int, default = 8, help="Number of parallel threads that pyABC will use")
    parser.add_argument("--n_pop", type=int, default=150, help="Population size, i.e. number of accepted samples (particles) to move on to new generation")
    parser.add_argument("--max_gen", type=int, default=10, help="Maximum number of generations to consider")
    parser.add_argument("--min_eps", type=float, default=0.01, help="Minimum acceptance threshold before ending the ABC Sampling prematurely")
    parser.add_argument("-G", "--genes_path", type=str, default="data/hg19_1000", help="Path to the genes list to be used.")
    parser.add_argument("-D", "--data_path", type=str, default="data/pcawg_filtered_95_pc.tsv", help="Path to the cancer cohort you want to match the fitness of.")
    parser.add_argument("--all_chromosomes", action="store_true", help="Flag to run SimChA with all chromosomes.")
    args = parser.parse_args()

    # Build the program once at the beginning
    subprocess.run(["dotnet build SimChA"], shell = True)
    
    # Initialize the input required for the fitness sampling.
    # These inputs are run-agnostic, so we just set it up once at the beginning
    cohort_path = args.data_path
    genes_path = args.genes_path
    inputs_path = "fitness_in"
    check_inputs(inputs_path, genes_path, cohort_path)

    bootstrap_path = os.path.join(inputs_path, "clones.tsv")
    binned_path    = os.path.join(inputs_path, "binned_CNs.tsv")

    # We use a wrapper function so that we can run SimChA with the relevant inputs and we only need to change them here
    def model_wrapper(params):
        return {"distance": run_simcha(params, genes_path, cohort_path, bootstrap_path, binned_path, args.repeats, args.all_chromosomes)}


    pwd = os.getcwd()
    out_dir = args.name
    subprocess.run([f"mkdir -p {out_dir}"], shell=True)

    # The fitness parameters (a, b, and c) are Dirichlet-distributed, so we keep c fixed
    # and normalize to 1 in the function "update_params_file"
    prior = Distribution(alpha=RV("uniform", 0, 1), beta=RV("uniform", 0, 1), w_strength = RV("uniform", 0.01, 25))
    # SimChA calculates the distance between simulated and observation, so we don't need an observed distance
    observed_data = {"distance": 0.0}
    sampler = sampler.MulticoreEvalParallelSampler(n_procs=args.n_procs)
    abc = ABCSMC(model_wrapper, prior, distance_function=distance, population_size = args.n_pop, sampler = sampler)
    # ABC-SMC output is a SQL database
    db_path = f"{out_dir}/parameters.db"
    abc.new("sqlite:///"+db_path, observed_data)

    # Main part of the ABC Sampling
    history = abc.run(minimum_epsilon = args.min_eps, max_nr_populations = args.max_gen)

    # Visualization of the posteriors of the fitness parameters
    df, w = history.get_distribution(m=0, t=history.max_t)
    arr_ax = plot_kde_matrix(df, w)
    fig = arr_ax[0,0].get_figure()
    fig.savefig(f"{out_dir}/kde_matrix.png", dpi=300, bbox_inches="tight")

    """fig, ax = plt.subplots()
    for t in range(history.max_t + 1):
        df, w = history.get_distribution(t=t)
        visualization.plot_kde_1d(
            df,
            w,
            xmin=30,
            xmax=150,
            x="n_events",
            xname=r"Event Count",
            ax=ax,
            label=f"PDF t={t}",
        )
    ax.legend()
    plt.savefig(f"{out_dir}/posterior_generations.png")"""


