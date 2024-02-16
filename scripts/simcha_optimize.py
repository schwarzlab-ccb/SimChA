#!/usr/bin/env python3
import subprocess
import numpy as np
import os
import argparse
import platform
from datetime import datetime
import uuid
import json

def generate_cohort(path, param_file, genes_path, out_path, repeats, all_chromomsomes):
    cmd = f"dotnet run --no-build --project {path} -C {param_file} -R {repeats} -O {out_path} -D {genes_path}"
    if not all_chromosomes:
        cmd += " --autosomes-only"
    subprocess.run([cmd], shell=True)
    return
   
def run_simcha(path, param_file, genes_path, cohort_path, repeats, all_chromosomes, out_dir):
    # Run command for SimChA
    cmd = f"dotnet run --no-build --project {path} -C {param_file} -R {repeats} -O {out_dir} --optimization -D {genes_path} -P {cohort_path}"
    if not all_chromosomes:
        cmd += " --autosomes-only"
    subprocess.run([cmd], shell=True)
   
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="pyABC program to fit parameters in SimChA ")
    parser.add_argument("-P", "--path", type=str, default="SimChA", help="Path to SimChA")
    parser.add_argument('-N', "--name", type=str, default="optimization_results", help="Name for output directory to put SQL database produced by pyABC and the posterior plot produced")
    parser.add_argument("-R", "--repeats", type=int, default=2000, help="Number of samples/repeats to generate for each SimChA run")
    parser.add_argument('--n_procs', type=int, default = 8, help="Number of parallel threads that pyABC will use")
    parser.add_argument("--n_pop", type=int, default=150, help="Population size, i.e. number of accepted samples (particles) to move on to new generation")
    parser.add_argument("-G", "--genes_path", type=str, default="data/hg19_1000", help="Path to the genes list to be used.")
    parser.add_argument("-D", "--data_path", type=str, default="data/pcawg_filtered_95_pc.tsv", help="Path to the cancer cohort you want to match the fitness of.")
    parser.add_argument("-T", "--test", action="store_true", help="Flag to use SimChA-generated data as the ground-truth rather than real data samples.")
    parser.add_argument("-C", "--config_file", type=str, default="simple_params.json", help="Parameter file used by SimChA.")
    parser.add_argument("--all_chromosomes", action="store_true", help="Flag to run SimChA with all chromosomes.")
    args = parser.parse_args()

    genes_path = args.genes_path
    # Build the program once for the corresponding thread, in case any changes have been made to SimChA
    subprocess.run([f"dotnet build {args.path}"], shell = True)

    # If we use SimChA-generated data as the ground-truth, we first have to generate that data
    param_file = args.config_file
    if args.test:
        cohort_dir_path = "out/ground_truth"
        # Generate the simulated data
        generate_cohort(args.path, param_file, genes_path, cohort_dir_path, args.repeats, args.all_chromosomes)
        cohort_path = "out/ground_truth/copynumbers.tsv"
    else:
        cohort_path = args.data_path
    
    # Create the output directory for the final plots as well as the SQL database
    out_dir = args.name
    subprocess.run([f"mkdir -p {out_dir}"], shell=True)

    run_simcha(args.path, param_file, genes_path, cohort_path, args.repeats, args.all_chromosomes, out_dir)
