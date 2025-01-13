# %%
import numpy as np
import json
from os.path import join
import subprocess
import uuid

# %%
# This python script is designed to let samples evolve for a given number of events
# and to vary the number of events so that we can observe the effects of time on
# the ploidy distributions

simcha_path = ".."
input_config_file = "../configs/config.json"
time_slices = [5, 10, 25, 50, 100, 200]
# %%
def main():
    for time in time_slices:
        run_simcha(input_config_file, time)

def run_simcha(config_file, time):
    with open(config_file, "r") as f:
        config = json.load(f)
    config["EvoParams"]["MaxTime"] = time
    config["EvoParams"]["WGDAccelerationFactor"] = 4
    id = str(uuid.uuid4())
    run_config = join(simcha_path, "scripts","temp", f"config_{id}.json")
    with open(run_config, "w") as f:
        json.dump(config, f, indent=4)
    cmd = f"dotnet run --no-build --project {simcha_path}/SimChA -C {run_config} -D {simcha_path}/data/hg19 -e -R 2000 -O {join(simcha_path, 'out', 'time_'+time)} --light"
    subprocess.run([cmd], shell=True, check=True)
    return id
# %%
if __name__ == "__main__":
    main()