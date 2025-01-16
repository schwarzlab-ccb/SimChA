# %%
import subprocess as sp
import os
import pandas as pd
import numpy as np
try:
    from tqdm import tqdm
except ImportError:
    print('could not import tqdm')
    tqdm = lambda x: x

# %%
# DOTNET = "dotnet"
DOTNET = "/usr/local/share/dotnet/dotnet"

base_cmd = [DOTNET, "run", "--project", "../SimChA", "--"]

# Test execution
out = sp.run(base_cmd + ["--help"], capture_output=True, text=True)
print(out.stderr)

# %%
signatures = [
    "0A", "0B", "1A", "1B", "2A", "2B", "2C", "3A", "3B", "3C", "3D", "3E",
    "1_25A_75B", "1_50A_50B", "1_75A_25B"
    ]
paths = { sig_name: f"./sigs/sig{sig_name}.json" for sig_name in signatures}
print(paths)
copynumbers_path = "./out/copynumbers.out"

# %%
def run_simulation(sig_list, reps, verbose=False):
    results = []
    mapping_key = {}
    for i in tqdm(range(reps)):
        returncode = 1
        id = i + 1
        select_sig = np.random.choice(sig_list)
        if verbose:
            print(f"{id}. Selected signature: {select_sig}")
        mapping_key[id] = select_sig
        event_count = 25 + np.random.binomial(100, 0.5)
        while returncode != 0:
            out = sp.run(base_cmd + ["-C", paths[select_sig], "-D", str(event_count)], capture_output=True, text=True)
            # if out.returncode != 0:
            #     raise Exception(f"Subprocess failed: {out.stderr}")
            returncode = out.returncode
        current_cn = pd.read_csv(copynumbers_path, sep="\t")
        current_cn["sample_id"] = id
        results.append(current_cn)
    return results, mapping_key

# %%
challenges = {
    "challenge0": ["0A", "0B"], "challenge0_1000": ["0A", "0B"],
    "challenge1": ["1A", "1B"], "challenge1_1000": ["1A", "1B"],
    "challenge2": ["2A", "2B", "2C"],
    "challenge3": ["3A", "3B", "3C", "3D", "3E"],
    "challenge1_mixed": ["1_25A_75B", "1_50A_50B", "1_75A_25B"]
    }

def run_challenge(challenge_name, reps, verbose=False):
    results, mapping_key = run_simulation(challenges[challenge_name], reps, verbose=verbose)
    df = pd.concat(results)
    df.to_csv(f"./{challenge_name}.tsv", index=False, sep="\t")
    # save mapping_key to a file
    with open(f"./{challenge_name}.key", "w") as file:
        for k, v in mapping_key.items():
            file.write(f"{k}, {v}\n")


# %%
run_challenge("challenge0", 100)

# %%
run_challenge("challenge1", 100)

# %%
run_challenge("challenge2", 100)

# %%
run_challenge("challenge3", 100)

# %%
run_challenge("challenge1_1000", 1000)

# %%
run_challenge("challenge0_1000", 1000)

# %%
run_challenge("challenge1_mixed_1000", 1000)

# %%
run_challenge("challenge1_mixed", 100)

# %%
for challenge in ['challenge0', 'challenge0_1000', 'challenge1', 'challenge1_1000']:

    data = pd.read_csv(f"./{challenge}.tsv", sep="\t")
    keys = pd.read_csv(f"./{challenge}.key", sep=",", header=None, names=["sample_id", "signature"])


    for key in keys['signature'].unique():
        data.loc[data['sample_id'].isin(keys.loc[keys['signature']==key, 'sample_id'].values)].to_csv(
            f"challenges_per_sig/{challenge}_{key}.tsv", index=False, sep="\t")


