#!/bin/bash
# Array of ThetaFitness values
theta_array=(5 10 20)
save_files(){
  outdir="out/"
  finaldir="distribution/"
  suffix="_theta_$1"
  
  events_base_name="events"
  mv "${outdir}${events_base_name}".tsv "${finaldir}${events_base_name}${suffix}".tsv

  samples_base_name="samples"
  mv "${outdir}${samples_base_name}".tsv "${finaldir}${samples_base_name}${suffix}".tsv

  copynumbers_base_name="copynumbers"
  mv "${outdir}${copynumbers_base_name}".tsv "${finaldir}${copynumbers_base_name}${suffix}".tsv

  clones_base_name="clones"
  mv "${outdir}${clones_base_name}".tsv "${finaldir}${clones_base_name}${suffix}".tsv

  karyotypes_base_name="karyotypes"
  mv "${outdir}${karyotypes_base_name}".tsv "${finaldir}${karyotypes_base_name}${suffix}".tsv

  sim_params="sim_params" 
  mv "${outdir}${sim_params}".json "${finaldir}${sim_params}${suffix}".json

}

for theta in "${theta_array[@]}"
do 
  # Modify the theta fitness parameter in the config file, using the template_params.json file
  # sed takes as input 1: where is the replacement, 2: old text, 3: new text
  sed_str='s/@@THETA@@/00/'
  cat template_params.json | sed -e "${sed_str/00/${theta}}"> my_params.json
  echo ${theta}
  # Do the simulation
  dotnet run --project SimChA -C my_params.json -M -R 2000
  # Save the end results
  save_files $theta
done
