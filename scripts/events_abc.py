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

    configs["EventCount"] = int(params["n_events"])
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

    # Build the program once for the corresponding thread
    subprocess.run(["dotnet build"], shell=True)
    
    def model_wrapper(params):
        return {"distance": run_simcha(params, genes_path, cohort_path, args.repeats)}


    pwd = os.getcwd()
    out_dir = "events_abc_results"
    subprocess.run([f"mkdir -p {out_dir}"], shell=True)

    # Uniform prior distributions for the various different properties of the simple events
    # We can also remove the number of events if we want
    # Event count is an integer so it has to be handled slightly differently
    event_count_values = np.arange(30,151)
    event_count_prob   = [1/len(event_count_values)]*len(event_count_values)

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


    prior = Distribution(n_events    = RV("rv_discrete", values=(event_count_values, event_count_prob)),
                         **{key: RV("uniform", a, b-a) for key, (a, b) in limits.items()})
    """w_chrom_del = RV("uniform", 1, 10),
    w_chrom_dup = RV("uniform", 1, 10), 
    w_int_dup   = RV("uniform", 1, 100),
    w_int_del   = RV("uniform", 1, 100),
    w_inv_dup   = RV("uniform", 1, 100),
    w_int_inv   = RV("uniform", 1, 100),
    w_bfb       = RV("uniform", 1, 10),
    w_tail_del  = RV("uniform", 1, 10),
    w_wgd       = RV("uniform", 0.01, 0.1),
    w_transloc  = RV("uniform", 1, 25),
    # Lengths are in units of 100kb
    l_int_dup   = RV("uniform", 1, 50),
    l_int_del   = RV("uniform", 1, 50),
    l_int_inv   = RV("uniform", 1, 50),
    l_inv_dup   = RV("uniform", 1, 50),
    l_transloc  = RV("uniform", 1, 50))"""
    event_limit = dict(n_events = (30, 150))
    limits.update(event_limit)

    transition = pyabc.AggregatedTransition(
	mapping={
		'n_events': pyabc.DiscreteJumpTransition(domain=event_count_values, p_stay=0.7),
		'w_chrom_del' : pyabc.MultivariateNormalTransition(),
		'w_chrom_dup' : pyabc.MultivariateNormalTransition(),
		'w_int_dup'   : pyabc.MultivariateNormalTransition(),
		'w_int_del'   : pyabc.MultivariateNormalTransition(),
		'w_inv_dup'   : pyabc.MultivariateNormalTransition(),
		'w_int_inv'   : pyabc.MultivariateNormalTransition(),
		'w_bfb'       : pyabc.MultivariateNormalTransition(),
		'w_tail_del'  : pyabc.MultivariateNormalTransition(),
		'w_wgd'       : pyabc.MultivariateNormalTransition(),
		'w_transloc'  : pyabc.MultivariateNormalTransition(),
		'l_int_del'   : pyabc.MultivariateNormalTransition(),
		'l_int_dup'   : pyabc.MultivariateNormalTransition(),
		'l_int_inv'   : pyabc.MultivariateNormalTransition(),
		'l_inv_dup'   : pyabc.MultivariateNormalTransition(),
		'l_transloc'  : pyabc.MultivariateNormalTransition()
		})	
	
    # SimChA calculates the Euclidean-summed Wasserstein distance, so we don't need an observed distance
    observed_data = {"distance": 0.0}
    sampler = sampler.MulticoreEvalParallelSampler(n_procs=16)
    abc = ABCSMC(model_wrapper, prior, distance_function=distance, transitions=transition, population_size = 150, sampler = sampler)
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
    plt.savefig(f"{out_dir}/posterior_event_count.png")
    df, w = history.get_distributions(m=0)
    ax = plot_kde_matrix(df, w, limits=limits)
    fig = ax.get_figure()
    fig.savefig(f"{out_dir}/kde_matrix.png", dpi=300, bbox_inches="tight")

