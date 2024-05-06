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
    if not all_chromomsomes:
        cmd += " --autosomes-only"
    subprocess.run([cmd], shell=True)
    return
   
def run_events_optimization(path, param_file, genes_path, cohort_path, repeats, all_chromosomes, out_dir, target_params):
    # Run command for SimChA
    cmd = f"dotnet run --no-build --project {path} -C {param_file} -R {repeats} -O {out_dir} --optimization -D {genes_path} -P {cohort_path} --target-params {target_params}"
    if not all_chromosomes:
        cmd += " --autosomes-only"
    subprocess.run([cmd], shell=True)

def run_fitness_optimization(path, param_file, genes_path, cohort_path, repeats, all_chromosomes, out_dir, binned_path, bootstrap_path):
    # Run command for SimChA
    cmd = f"dotnet run --no-build --project {path} -C {param_file} -R {repeats} -O {out_dir} --optimization fitness -D {genes_path} -P {cohort_path} --binned-samples {binned_path} -B {bootstrap_path} -M"
    if not all_chromosomes:
        cmd += " --autosomes-only"
    subprocess.run([cmd], shell=True)

def check_inputs(in_path, genes_path, cohort_path):

    clones_file = "clones.tsv"
    binned_file = "binned_CNs.tsv"

    if not os.path.isdir(in_path):
        subprocess.run([f"mkdir -p in_path"], shell = True)
    # bootstrap file needed for sampling the fitnesses of the mcmc produces clones
    if not os.path.isfile(os.path.join(in_path, clones_file)):
        cmd = f"dotnet run --no-build --project SimChA -P {cohort_path} -O {in_path} -D {genes_path}"
        subprocess.run([cmd], shell = True)
    # To avoid having to generate the binned copy-number profiles for each run of SimChA, we do it once here if the file does not exist.
    if not os.path.isfile(os.path.join(in_path, binned_file)):
        cmd = f"dotnet run --no-build --project SimChA -P {cohort_path} -O {in_path} --bin-samples"
        subprocess.run([cmd], shell = True)

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
    parser.add_argument("--fitness", action="store_true", help="Flag to run the fitness optimization.")
    parser.add_argument("--target_params", type=str, default="target_params.json", help="Target parameters config file")
    args = parser.parse_args()

    genes_path = args.genes_path
    # Build the program once for the corresponding thread, in case any changes have been made to SimChA
    subprocess.run([f"dotnet build {args.path}"], shell = True)

    # If we use SimChA-generated data as the ground-truth, we first have to generate that data
    param_file = args.config_file
    
    # Create the output directory for the final plots as well as the SQL database
    out_dir = args.name
    subprocess.run([f"mkdir -p {out_dir}"], shell=True)

    if args.fitness:
        cohort_path = args.data_path
        inputs_path = "fitness_in"
        check_inputs(inputs_path, genes_path, cohort_path)
        bootstrap_path = os.path.join(inputs_path, "clones.tsv")
        binned_path    = os.path.join(inputs_path, "binned_CNs.tsv")
        run_fitness_optimization(args.path, param_file, genes_path, cohort_path, args.repeats, args.all_chromosomes, out_dir, binned_path, bootstrap_path)
    else:
        if args.test:
            cohort_dir_path = "out/ground_truth"
            # Generate the simulated data
            generate_cohort(args.path, args.target_params, genes_path, cohort_dir_path, args.repeats, args.all_chromosomes)
            cohort_path = "out/ground_truth/copynumbers.tsv"
        else:
            cohort_path = args.data_path
        run_events_optimization(args.path, param_file, genes_path, cohort_path, args.repeats, args.all_chromosomes, out_dir, args.target_params)
