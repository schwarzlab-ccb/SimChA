#!/usr/bin/env python3
import subprocess
import os
import argparse
import platform
import datetime as dt
import json
from pyabc import ABCSMC, RV, Distribution, settings
from pyabc.populationstrategy import ConstantPopulationSize

settings.set_figure_params('pyabc')

def update_params_file(params):
    # Create the temporary parameter file
    foldername = f"{int(dt.datetime.now().timestamp())}"
    path = f"{pwd}/temp/{foldername}"
    subprocess.run([f"mkdir -p {path}"], shell = True)
    subprocess.run([f"cp simple_params.json {path}"], shell = True)

    file_path = f"{path}/simple_params.json"

    # Update the parameter file
    with open(file_path, 'r', encoding='utf-8') as json_file:
        configs = json.load(json_file)

    configs["EventCount"] = params["n_events"]

    """configs["Signatures"]["CNs"]["Events"]["ChromDeletion"]["Prob"]        = params["w_chrom_del"]
    configs["Signatures"]["CNs"]["Events"]["ChromDuplication"]["Prob"]     = params["w_chrom_dup"]
    configs["Signatures"]["CNs"]["Events"]["InternalDuplication"]["Prob"]  = params["w_int_dup"]
    configs["Signatures"]["CNs"]["Events"]["InternalDeletion"]["Prob"]     = params["w_int_del"]
    configs["Signatures"]["CNs"]["Events"]["InternalInversion"]["Prob"]    = params["w_int_inv"]
    configs["Signatures"]["CNs"]["Events"]["InvertedDuplication"]["Prob"]  = params["w_inv_dup"]
    configs["Signatures"]["CNs"]["Events"]["BreakageFusionBridge"]["Prob"] = params["w_bfb"]
    configs["Signatures"]["CNs"]["Events"]["TailDeletion"]["Prob"]         = params["w_tail_del"]
    configs["Signatures"]["CNs"]["Events"]["WholeGenomeDoubling"]["Prob"]  = params["w_wgd"]
    
    configs["Signatures"]["CNs"]["Events"]["InternalInversion"]["Size"]   = params["l_int_inv"]*100_000
    configs["Signatures"]["CNs"]["Events"]["InternalDeletion"]["Size"]    = params["l_int_del"]*100_000
    configs["Signatures"]["CNs"]["Events"]["InternalDuplication"]["Size"] = params["l_int_dup"]*100_000
    configs["Signatures"]["CNs"]["Events"]["InvertedDuplication"]["Size"] = params["l_inv_dup"]*100_000"""
    with open(file_path, 'w', encoding="utf-8") as json_file:
        json.dump(configs, json_file)
    # Return the path to the config file
    return file_path


def run_simcha(params):
    param_file_path = update_params_file(params)

    cmd = ["dotnet", "run", "--project",  "SimChA", "-C", f"{param_file_path}" "-R", "2000", "-O", f"abc_results", "--optimization", "-D", "data/hg19_1000", "-P", "pcawg_filtered_95_pc.tsv"]
    output = subprocess.check_output(cmd, universal_newlines=True, shell=True)
    # SimChA produces as its output the Euclidean sum of Wasserstein distances for each of the 
    # characteristic features of cancer genomes, printing the double to the command 
    last_line = output.strip().split("\n")[-1]
    return float(last_line.split(":")[1].strip())

def model(params):
    return {"distance": run_simcha(params)}

def distance(x,y):
    return abs(x["distance"] - y["distance"])
    

if __name__ == "__main__":
    pwd = os.getcwd()
    out_dir = "abc_results"
    subprocess.run([f"mkdir -p {out_dir}"], shell=True)

    # Uniform prior distributions for the various different properties of the simple events
    # We can also remove the number of events if we want
    prior = Distribution(n_events    = RV("uniform", 30, 150))
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

    # SimChA calculates the Euclidean-summed Wasserstein distance, so we don't need an observed distance
    observed_data = {"distance": 0.0}

    abc = ABCSMC(model, prior, distance_function=distance, population_size = 5)
    # ABC-SMC output is a SQL database
    db_path = f"{out_dir}/test.db"
    abc.new("sqlite:///"+db_path, observed_data)

    history = abc.run(minimum_epsilon = 0.1, max_nr_populations = 1)

    fig, ax = plt.subplots()
    for t in range(history.max_t + 1):
        df, w = history.get_distribution(t=t)
        pyabc.visualization.plot_kde_1d(
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


