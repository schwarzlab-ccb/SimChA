#!/usr/bin/env nextflow
nextflow.enable.dsl=2
simcha_path = workflow.launchDir + "/SimChA"
params.simcha_params_file = "/projects/ag-schwarzr/project-simcha/simcha/configs/main_config.json"

process SimChA {
	publishDir "${workflow.launchDir}/results_event_cost/${e_cost}_${alpha}", mode: 'move'
	
	input:
	tuple val(e_cost), val(alpha)
	
	output:
	path("*")	

	script:

	"""
	cp ${params.simcha_params_file} config.json
	~/.conda/envs/simcha/bin/python -c "
	import json
	import numpy as np
	def calc_pwgd(alpha):
		r = 0.35 * np.exp(alpha) / 89.3
		return r/(1-r)
	with open('config.json', 'r') as f:
		config = json.load(f)
	config['Fitness']['Stress'] = ${alpha}
	config['Fitness']['TsgOg'] = 0
	config['Fitness']['Essentiality'] = 0
	config['Signatures']['CNVs']['Events'][-1]['Prob'] = calc_pwgd(${alpha})
	config['EvoParams']['EventCost'] = ${e_cost}	
	config['EvoParams']['TetraploidStart'] = False
	with open('config.json', 'w') as f:
		json.dump(config, f, indent=4)
	"
	~/.conda/envs/simcha/lib/dotnet/dotnet run --no-build --project ${simcha_path} -- -C config.json -D ${workflow.launchDir}/data/hg19 -e -R 10000 -O "." --light
        """
}

workflow {
	def max_vals_per_param = 50
	def e_cost = Channel.from(params.e_cost).take(100)
	def alpha = Channel.from(params.alpha).take(max_vals_per_param)
	def product = e_cost.combine(alpha)
	SimChA(product)
}
