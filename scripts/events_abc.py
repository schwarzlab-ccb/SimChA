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

    configs["EventCount"] = int(params["n_events"])

    """configs["Signatures"]["CNs"]["Events"]["ChromDeletion"]["Prob"]        = params["w_chrom_del"]
    configs["Signatures"]["CNs"]["Events"]["ChromDuplication"]["Prob"]     = params["w_chrom_dup"]
    configs["Signatures"]["CNs"]["Events"]["InternalDuplication"]["Prob"]  = params["w_int_dup"]
    configs["Signatures"]["CNs"]["Events"]["InternalDeletion"]["Prob"]     = params["w_int_del"]
    configs["Signatures"]["CNs"]["Events"]["InternalInversion"]["Prob"]    = params["w_int_inv"]
    configs["Signatures"]["CNs"]["Events"]["InvertedDuplication"]["Prob"]  = params["w_inv_dup"]
    configs["Signatures"]["CNs"]["Events"]["BreakageFusionBridge"]["Prob"] = params["w_bfb"]
    configs["Signatures"]["CNs"]["Events"]["TailDeletion"]["Prob"]         = params["w_tail_del"]
    configs["Signatures"]["CNs"]["Events"]["WholeGenomeDoubling"]["Prob"]  = params["w_wgd"]
    # The length-scales have to be rounded to integer status 
    configs["Signatures"]["CNs"]["Events"]["InternalInversion"]["Size"]   = int(params["l_int_inv"]*100_000)
    configs["Signatures"]["CNs"]["Events"]["InternalDeletion"]["Size"]    = int(params["l_int_del"]*100_000)
    configs["Signatures"]["CNs"]["Events"]["InternalDuplication"]["Size"] = int(params["l_int_dup"]*100_000)
    configs["Signatures"]["CNs"]["Events"]["InvertedDuplication"]["Size"] = int(params["l_inv_dup"]*100_000)"""
    with open(file_path, 'w', encoding="utf-8") as json_file:
        json.dump(configs, json_file)
    # Return the path to the config file
    return path


def run_simcha(params, genes_path, cohort_path, repeats):
    param_file_path = update_params_file(params)

    cmd = f"dotnet run --no-build --project SimChA -C {param_file_path}/simple_params.json -R {repeats} -O {param_file_path}/out --optimization events -D {genes_path} -P {cohort_path}"
    output = subprocess.check_output([cmd], universal_newlines=True, shell=True)
    # SimChA produces as its output the Euclidean sum of Wasserstein distances for each of the 
    # characteristic features of cancer genomes, printing the double to the command 
    last_line = output.strip().split("\n")[-1]

    # Delete the temporary folder and files
    subprocess.run([f"rm -rf {param_file_path}"], shell=True)
    # Return the distance SimChA calculated
    return float(last_line.split(":")[1].strip()

def distance(x,y):
    return abs(x["distance"] - y["distance"])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="pyABC program to fit parameters in SimChA ")
    #parser.add_argument('-R', "--repeats", type=int, default=500, help="Number of SimChA simulated samples to generate for each pyABC sample")
    parser.add_argument('-N', "--name",type=str, default="", help="Name for output directory to put SQL database produced by pyABC and the posterior plot produced")
    parser.add_argument("-r", "--repeats", type=int, default=2000, help="Number of samples to generate for each SimChA run")
    args = parser.parse_args()


    cohort_path = "data/pcawg_filtered_95_pc.tsv"
    genes_path = "data/hg19_1000"

    def model_wrapper(params):
        return {"distance": run_simcha(params, genes_path, cohort_path, args.repeats)}


    pwd = os.getcwd()
    out_dir = "events_abc_results"
    subprocess.run([f"mkdir -p {out_dir}"], shell=True)

    # Uniform prior distributions for the various different properties of the simple events
    # We can also remove the number of events if we want
    event_count_values = np.arange(30,151)
    event_count_prob   = [1/len(event_count_values)]*len(event_count_values)


    prior = Distribution(n_events    = RV("rv_discrete", values=(event_count_values, event_count_prob)))
    """w_chrom_del = RV("uniform", 1, 10), 
            w_chrom_dup = RV("uniform", 1, 10), 
            w_int_dup   = RV("uniform", 1, 100),
            w_int_del   = RV("uniform", 1, 100),
            w_inv_dup   = RV("uniform", 1, 100),
            w_int_inv   = RV("uniform", 1, 100),
            w_bfb       = RV("uniform", 1, 10),
            w_tail_del  = RV("uniform", 1, 10),
            w_wgd       = RV("uniform", 0.01, 0.1),
            # Lengths are in units of 100kb
            l_int_dup   = RV("uniform", 1, 50),
            l_int_del   = RV("uniform", 1, 50),
            l_int_inv   = RV("uniform", 1, 50),
            l_inv_dup   = RV("uniform", 1, 50)
            )"""

    transition = pyabc.AggregatedTransition(
	mapping={
		'n_events': pyabc.DiscreteJumpTransition(domain=event_count_values, p_stay=0.7)
		}
	)
    # Transitions for the other parameters:
    """
	'w_chrom_del' : pyabc.MultivariateNormalTransition(),
	'w_chrom_dup' : pyabc.MultivariateNormalTransition(),
	'w_int_dup'   : pyabc.MultivariateNormalTransition(),
	'w_int_del'   : pyabc.MultivariateNormalTransition(),
	'w_inv_dup'   : pyabc.MultivariateNormalTransition(),
	'w_int_inv'   : pyabc.MultivariateNormalTransition(),
	'w_bfb'       : pyabc.MultivariateNormalTransition(),
	'w_tail_del'  : pyabc.MultivariateNormalTransition(),
	'w_wgd'       : pyabc.MultivariateNormalTransition(),
	'l_int_del'   : pyabc.MultivariateNormalTransition(),
	'l_int_dup'   : pyabc.MultivariateNormalTransition(),
	'l_int_inv'   : pyabc.MultivariateNormalTransition(),
	'l_inv_dup'   : pyabc.MultivariateNormalTransition()
	
	
    """

    # SimChA calculates the Euclidean-summed Wasserstein distance, so we don't need an observed distance
    observed_data = {"distance": 0.0}
    sampler = sampler.MulticoreEvalParallelSampler(n_procs=16)
    abc = ABCSMC(model_wrapper, prior, distance_function=distance, transitions=transition, population_size = 100, sampler = sampler)
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


