import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from os.path import join
import BISCUT_preprocessing as pre
import os

# Autosomes
chrs = [f"chr{i}" for i in range(1, 23)]
chr_arms = {i: ["p", "q"] for i in chrs}
# ignore acrocentric chromosomes
chr_arms["chr13"] = ["q"]
chr_arms["chr14"] = ["q"]
chr_arms["chr15"] = ["q"]
chr_arms["chr21"] = ["q"]
chr_arms["chr22"] = ["q"]

telcent_threshold = 10e-3

def preprocess_data(name):
    data = pd.read_csv(f"../{name}/copynumbers.tsv", sep="\t")
    data["Num_Probes"] = 10
    data["Segment_Mean"] = (data["cn_a"] + data["cn_b"])/2.0
    data.loc[data["Segment_Mean"] < 0.00001, 'Segment_Mean'] = 0.1
    data["Segment_Mean"] = np.log2(data["Segment_Mean"])
    data = data.drop(["cn_a", "cn_b", "n_snvs"], axis=1)
    data["chrom"] = data["chrom"].str[3:]
    data = data.rename(columns={"chrom": "Chromosome", "sample_id" : "Sample", "start":"Start", "end":"End"})
    data = data[~data["Chromosome"].isin(["X","Y"])]
    data.to_csv(f"docs/{name}.tsv", sep="\t", index=False)
    pre.take_care_arms(name, "2023_07_27")

def plot_partial_scna(output, preprocess):
    names = ["pcawg", "simple_no_fit", "complex_no_fit","mcmc_simple", "mcmc_complex"]
    if preprocess:
        for name in names:
            preprocess_data(name)
    print(names)
    types = ["amp_cent", "amp_tel", "del_cent", "del_tel"]

    for i, chr_i in enumerate(chrs):
        n_arms = len(chr_arms[chr_i])
        print(f"arm numbers: {n_arms}")
        fig, ax = plt.subplots(n_arms, 4, figsize=(16*len(types),9*n_arms))
        # acrocentric chromosomes:
        for k in range(n_arms):
            if n_arms == 1:
                arm_name = ""
            else:
                arm_name = "p" if k == 0 else "q"
            for j, scna_type in enumerate(types):
                current_type = f"{chr_i[3:]}{arm_name}_{scna_type}"
                hist = []
                for name in names:
                    dir = f"breakpoint_files_2023_07_27/{name}"
                    arm_results = pd.read_csv(f"{dir}/{name}_{current_type}.txt", sep="\t")
                    # Filter the arm-results:
                    arm_results = arm_results[arm_results["percent"] >= telcent_threshold]
                    arm_results = arm_results[arm_results["percent"] <= (1.0-telcent_threshold)]
                    hist.append(arm_results["percent"])
                if n_arms == 1:
                    ax[j].violinplot(hist)
                    ax[j].set_title(f"{current_type}")
                    ax[j].set_xticks(range(1,6), ["pcawg", "simple", "complex", "mcmc simp", "mcmc comp"])
                else:
                    ax[k,j].violinplot(hist)
                    ax[k,j].set_title(f"{current_type}")
                    ax[k,j].set_xticks(range(1,6),["pcawg", "simple", "complex", "mcmc simp", "mcmc comp"])
        fig.savefig(join(output, f"{chr_i}.png"), dpi=150, bbox_inches='tight')
        plt.close()
    #filenames = os.listdir()
    #for f in filenames:
    
    #hist = get_partial_scna(data, centromeres)

    #ax.violinplot(hist)
    #ax.set_ylim(0,1)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot the partial SCNA violin plots for a number of datasets')
    parser.add_argument('-i', "--input", type=str, help='The folder to input dataset to plot')
    parser.add_argument("-p", "--preprocess",action="store_true", default=False, help="Preprocess the data")
    parser.add_argument('-o', "--output", type=str, help='The folder to output the plots to')
    args = parser.parse_args()

    plot_partial_scna(args.output, args.preprocess)