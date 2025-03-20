#!/usr/bin/env nextflow
nextflow.enable.dsl=2
simcha_path = workflow.launchDir + "/SimChA"
params.simcha_params_file = "/projects/ag-schwarzr/project-simcha/simcha/configs/C_ISMB_fit.json"

process SimChA {
	publishDir "${workflow.launchDir}/results_ISMB_ess_scan/${ess}", mode: 'move'
	
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
	config['Fitness']['TsgOg'] = 0
	config['Fitness']['Essentiality'] = ${ess}
	config['Signatures']['CNVs']['Events'][-1]['Prob'] = 0
	config['EvoParams']['PWGD'] = 0
	config['Fitness']['Delta'] = 0
	with open('config.json', 'w') as f:
		json.dump(config, f, indent=4)
	"
	~/.conda/envs/simcha/lib/dotnet/dotnet run --no-build --project ${simcha_path} -- -C config.json -D ${workflow.launchDir}/data/hg19 -e -R 10000 -O "." --light
        """
}

workflow {
	def max_vals_per_param = 1000
	def ess = Channel.from(params.ess_1000).take(max_vals_per_param)
	SimChA(ess)
}
