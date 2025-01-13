#!/usr/bin/env nextflow
nextflow.enable.dsl=2
simcha_path = workflow.launchDir + "/SimChA"
params.simcha_params_file = "/projects/ag-schwarzr/project-simcha/simcha/configs/all_events_WGD.json"

process SimChA {
	publishDir "${workflow.launchDir}/results_all_events_WGD/${ess}", mode: 'move'
	
	input:
	val ess
	
	output:
	path("*")
	
	script:

	"""
	cp ${params.simcha_params_file} config.json
	~/.conda/envs/simcha/bin/python -c "
	import json
	with open('config.json', 'r') as f:
		config = json.load(f)
	config['Fitness']['Essentiality'] = ${ess}
	with open('config.json', 'w') as f:
		json.dump(config, f, indent=4)
        "
	~/.conda/envs/simcha/lib/dotnet/dotnet run --no-build --project ${simcha_path} -- -C config.json -D ${workflow.launchDir}/data/hg19 -e -R 10000 -O "." --light
        """
}

workflow {
	def max_vals_per_param = 500
	def ess = Channel.from(params.ess).take(max_vals_per_param)
	SimChA(ess)
}
