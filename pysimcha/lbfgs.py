#!/usr/bin/env python3
import numpy as np
import pandas as pd
import subprocess
import argparse
import platform
import json
import os
from scipy.optimize import minimize
from scipy.stats import gaussian_kde
#from utils import get_seg_len_means

def get_seg_len_means(data, include_cn_normal = False, include_loh = False, include_sex_chromosomes = False):
    sample_ids = pd.Series({c: data[c].unique() for c in data})["sample_id"]
    all_means = []
    for _, id in enumerate(sample_ids):
        df = data[data["sample_id"] == id]
        # Remove sex chromosomes
        if not include_sex_chromosomes:
            df = df[~df["chrom"].isin(["chrX", "chrY"])]
        # Remove cn normal and LOH segments if the flags are turned on
        if not include_cn_normal:
            df = df[~((df["cn_a"] == 1) & (df["cn_b"] == 1))]
        elif not include_loh:
            df = df[~((df["cn_a"] != 1) & (df["cn_a"] + df["cn_b"] == 2))]
        
        all_means.append((df['end'] - df['start']).mean() if not df.empty else 0)
    return all_means

def update_params(params, out_dir):
    file_path = os.path.join(out_dir, param_file)
    with open(file_path, 'r', encoding="utf-8") as json_file:
        configs = json.load(json_file)

    # ChromDeletion
    configs["Signatures"]["CNs"]["Events"][0]["Prob"] = params[0]
    # ChromDuplication
    configs["Signatures"]["CNs"]["Events"][1]["Prob"] = params[1] 
    # InternalDuplication
    configs["Signatures"]["CNs"]["Events"][2]["Prob"] = params[2] 
    # InternalDeletion
    configs["Signatures"]["CNs"]["Events"][3]["Prob"] = params[3] 
    # BreakageFusionBridge    
    configs["Signatures"]["CNs"]["Events"][6]["Prob"] = params[4] 
    # TailDeletion
    configs["Signatures"]["CNs"]["Events"][7]["Prob"] = params[5]
    # Whole-Genome Doubling
    configs["Signatures"]["CNs"]["Events"][9]["Prob"] = 1.0 - np.sum(params)

    with open(file_path, "w", encoding="utf-8") as json_file:
        json.dump(configs, json_file)
    return file_path

def generate_cohort(param_file, genes_path, out_path, repeats, all_chromomsomes):
    cmd = f"dotnet run --no-build --project SimChA -C {param_file} -R {repeats} -O {out_path} -D {genes_path}"
    if not all_chromosomes:
        cmd += " --autosomes-only"
    subprocess.run([cmd], shell=True)
    return
   
def run_simcha(params, genes_path, cohort_path, repeats, all_chromosomes, out_dir):
    param_file = update_params(params, out_dir)
    # Run command for SimChA
    cmd = f"dotnet run --no-build --project SimChA -C {param_file} -R {repeats} -O {out_dir} -D {genes_path}"
    if not all_chromosomes:
        cmd += " --autosomes-only"
    subprocess.run([cmd], shell=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="pyABC program to fit parameters in SimChA ")
    parser.add_argument('-O', "--out_dir",type=str, default="events_abc_results", help="Name for output directory to put SQL database produced by pyABC and the posterior plot produced")
    parser.add_argument("-R", "--repeats", type=int, default=2000, help="Number of samples/repeats to generate for each SimChA run")
    parser.add_argument("-G", "--genes_path", type=str, default="data/hg19_100", help="Path to the genes list to be used.")
    parser.add_argument("-D", "--data_path", type=str, default="data/pcawg_filtered_95_pc.tsv", help="Path to the cancer cohort you want to match the fitness of.")
    parser.add_argument("-T", "--test", action="store_true", help="Flag to use SimChA-generated data as the ground-truth rather than real data samples.")
    parser.add_argument("-C", "--config_file", type=str, default="simple_params.json", help="Parameter file used by SimChA.")
    parser.add_argument("--all_chromosomes", action="store_true", help="Flag to run SimChA with all chromosomes.")
    args = parser.parse_args()

    # Build the program once for the corresponding thread, in case any changes have been made to SimChA
    subprocess.run(["dotnet build SimChA"], shell = True)

    # If we use SimChA-generated data as the ground-truth, we first have to generate that data
    param_file = args.config_file
    if args.test:
        cohort_dir_path = "out/ground_truth"
        # Generate the simulated data
        generate_cohort(param_file, args.genes_path, cohort_dir_path, args.repeats, args.all_chromosomes)
        cohort_path = "out/ground_truth/copynumbers.tsv"
    else:
        cohort_path = args.data_path

    # Create the output directory for the final plots as well as the SQL database
    subprocess.run([f"mkdir -p {args.out_dir}"], shell=True)
    subprocess.run([f"cp {args.config_file} {args.out_dir}"], shell=True)
    
    print("Getting Summary Feature distribution")
    # Generate the histogram used by the L-BFGS method
    df = pd.read_csv(cohort_path, sep="\t")
    data = get_seg_len_means(df)
    # Estimate PDF with KDE
    kde = gaussian_kde(data)

    # Define the negative log-likelihood function
    def neg_log_likelihood(params):
        print("new sample")
         # Ensure parameters don't violate constraints
        if any(param < 0 or param > 1 for param in params) or sum(params) > 1:
            print("error in params")
            print(params)
            return np.inf

        print(params)
        print(1.0 - np.sum(params))

        run_simcha(params, args.genes_path, cohort_path, args.repeats, args.all_chromosomes, args.out_dir)
        sim_df = pd.read_csv(os.path.join(args.out_dir, "copynumbers.tsv"), sep="\t")
        sim_data = get_seg_len_means(sim_df)
        # Change the segment length mean into units of MB
        sim_data = [l/1_000_000 for l in sim_data]
        # Calculate probabilities for each data point under the KDE
        probabilities = max(kde.evaluate(sim_data), 1e-10)
        print(probabilities[0])
        log_likelihood = np.sum(np.log(probabilities))
        return -log_likelihood

    # Define bounds to ensure non-negativity and sum-to-one
    def constraint(params):
        return 1 - sum(params)

    bounds = [(0, 1) for _ in range(6)]
    cons = {'type': 'eq', 'fun': constraint}

    # Initial params for simulated ground-truth
    initial_params = [10, 10, 50, 50, 10, 10]
    # Initial params for PCAWG
    # initial_params = [0.0058, 0.0054, 0.6755, 0.1983, 0.0697, 0.0453]#, 0.001]
    initial_params = [0.99*p/np.sum(initial_params) for p in initial_params]
    
    print("Starting optimization")
    # Optimize using a method that supports constraints
    result = minimize(neg_log_likelihood, initial_params, method='SLSQP', bounds=bounds, constraints=cons)
    wgd = 1.0 - np.sum(result.x)
    print(result.x)
    print(wgd)
    

