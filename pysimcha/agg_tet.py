import numpy as np
import pandas as pd
import cns
import os
from os.path import join as pjoin
import subprocess
import sys

# %%
def run_cns_aggregate(base_path, exp_path, segments_path):
    filepath = pjoin(base_path, exp_path, "copynumbers.tsv")
    out_path = pjoin(base_path, exp_path, "cns_3MB.tsv")
    if os.path.exists(out_path):
        print(f"File already exists: {out_path}")
        return 0
    print(f"Running cns aggregate for {filepath}")
    cmd = [
        'cns', 
        'aggregate',
        filepath,
        '--segments', segments_path,
        '--threads', "20",
        '--out', out_path,
        '--verbose'
    ]
        
    try:
        result = subprocess.run(
            cmd,
            check=False,
            capture_output=True,
            text=True
        )
        
        print(f"Exit code: {result.returncode}")
        
        if result.stderr:
            print("Errors:", result.stderr, file=sys.stderr)
        
        print("Result in:", out_path)
        return result.returncode
        
    except Exception as e:
        print(f"Failed to run command: {e}", file=sys.stderr)
        return 1

# %%
wgd_status = [0, 1]

#%%

curr_dir = pjoin("../", "results_ISMB_tsg_tetraploid_and_delta_scan")
segments = "segs_3MB.bed"

subdirs = [d for d in os.listdir(curr_dir) if os.path.isdir(pjoin(curr_dir, d))]
for subdir in subdirs:
    run_cns_aggregate(curr_dir, subdir, segments)
