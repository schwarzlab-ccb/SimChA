#!/usr/bin/env nextflow
nextflow.enable.dsl=2
simcha_path = workflow.launchDir + "/SimChA"
params.simcha_params_file = "/projects/ag-schwarzr/project-simcha/simcha/configs/C_ISMB_fit.json"

process SimChA {
	publishDir "${workflow.launchDir}/results_ISMB_tsg_tetraploid_and_delta_scan/${tsg}_${delta}", mode: 'move'
	
	input:
	tuple val(tsg), val(delta)
	
	output:
	path("*")
	
	script:
	"""
	# Check if the folder and the target file exist
	if [ -d "${workflow.launchDir}/results_ISMB_tsg_tetraploid_and_delta_scan/${tsg}_${delta}" ] && [ -f "${workflow.launchDir}/results_ISMB_tsg_tetraploid_and_delta_scan/${tsg}_${delta}/copynumbers.tsv" ]; then
		echo "Skipping process "
		touch "./skipped.txt"
        	exit 0
	fi
	cp ${params.simcha_params_file} config.json
	~/.conda/envs/simcha/bin/python -c "
	import json 
	with open('config.json', 'r') as f:
		config = json.load(f)
	config['Fitness']['Stress'] = 0
	config['Fitness']['TsgOg'] = ${tsg}
	config['Fitness']['Essentiality'] = 0
	config['Fitness']['Delta'] = ${delta}
	config['Signatures']['CNVs']['Events'][-1]['Prob'] = 0
	config['EvoParams']['PWGD'] = 0
	config['EvoParams']['TetraploidStart'] = True
	config['EvoParams']['EventCost'] = 1.6
	config['EvoParams']['MaxTries'] = 100
	with open('config.json', 'w') as f:
		json.dump(config, f, indent=4)
	"
	~/.conda/envs/simcha/lib/dotnet/dotnet run --no-build --project ${simcha_path} -- -C config.json -D ${workflow.launchDir}/data/hg19 -e -R 1000 -O "."
        """
}

workflow {
	def tsg = Channel.from(params.tsg).take(50)
	def delta = Channel.from(params.delta).take(50)
	def product = tsg.combine(delta)
	SimChA(product)
}
