#!/usr/bin/env python3
import subprocess
import numpy as np
import os
import argparse
import platform
import datetime as dt
import json
import matplotlib.pyplot as plt
from pyabc import ABCSMC, RV, Distribution, settings, visualization, sampler
from pyabc.populationstrategy import ConstantPopulationSize
import pyabc

settings.set_figure_params('pyabc')

def update_params_file(params):
    # Create the temporary parameter file
    
    foldername = f"{int(dt.datetime.now().timestamp())}_"+"_".join(str(p) for p in params.values())
    
    path = f"{pwd}/temp/{foldername}"
    subprocess.run([f"mkdir -p {path}"], shell = True)
    subprocess.run([f"cp simple_params.json {path}"], shell = True)
    # Create the output subdirectory
    subprocess.run([f"mkdir -p {path}/out"], shell=True)
    file_path = f"{path}/simple_params.json"

    # Update the parameter file
    with open(file_path, 'r', encoding='utf-8') as json_file:
        configs = json.load(json_file)
    
    configs["Fitness"]["Stress"], configs["Fitness"]["TsgOg"], configs["Fitness"]["Essentiality"] = params["abc"]

    with open(file_path, 'w', encoding="utf-8") as json_file:
        json.dump(configs, json_file)
    # Return the path to the config file
    return path


def run_simcha(params):
    param_file_path = update_params_file(params)

    cmd = f"dotnet run --no-build --project SimChA -C {param_file_path}/simple_params.json -R 250 -O {param_file_path}/out --optimization -D data/hg19_1000 -M -P pcawg_filtered_95_pc.tsv"
    output = subprocess.check_output([cmd], universal_newlines=True, shell=True)
    # SimChA produces as its output the Euclidean sum of Wasserstein distances for each of the 
    # characteristic features of cancer genomes, printing the double to the command 
    last_line = output.strip().split("\n")[-1]

    # Delete the temporary folder and files
    subprocess.run([f"rm -rf {param_file_path}"], shell=True)
    # Return the distance SimChA calculated
    return float(last_line.split(":")[1].strip())

def model(params):
    return {"distance": run_simcha(params)}

def distance(x,y):
    return abs(x["distance"] - y["distance"])
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="pyABC program to fit parameters in SimChA ")
    #parser.add_argument('-R', "--repeats", type=int, default=500, help="Number of SimChA simulated samples to generate for each pyABC sample")
    parser.add_argument('-N', "--name",type=str, default="", help="Name for output directory to put SQL database produced by pyABC and the posterior plot produced")
    parser.add_argument("-C", "--nCPUs", type=int, default=16, help="Number of cpus to use")

    args = parser.parse_args()


    pwd = os.getcwd()
    out_dir = "fitness_abc_results"#+args.name
    subprocess.run([f"mkdir -p {out_dir}"], shell=True)

    # Uniform prior distributions for the various different properties of the simple events
    # We can also remove the number of events if we want
    event_count_values = np.arange(30,151)
    event_count_prob   = [1/len(event_count_values)]*len(event_count_values)

    # Using symmetric concentration parameters as an example
    prior = Distribution(abc=RV("dirichlet", [1, 1, 1])

    # SimChA calculates the distance between simulated and observation, so we don't need an observed distance
    observed_data = {"distance": 0.0}
    sampler = sampler.MulticoreEvalParallelSampler(n_procs=args.nCPUs())
    abc = ABCSMC(model, prior, distance_function=distance, transitions=transition, population_size = 100, sampler = sampler)
    # ABC-SMC output is a SQL database
    db_path = f"{out_dir}/test.db"
    abc.new("sqlite:///"+db_path, observed_data)

    history = abc.run(minimum_epsilon = 0.01, max_nr_populations = 10)

    fig, ax = plt.subplots()
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
    plt.savefig(f"{out_dir}/posterior_generations.png")


