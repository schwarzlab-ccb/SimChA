#!/usr/bin/env nextflow
nextflow.enable.dsl=2
simcha_path = workflow.launchDir + "/SimChA"
params.simcha_params_file = "/projects/ag-schwarzr/project-simcha/simcha/configs/C_ISMB_fit.json"

process SimChA {
	publishDir "${workflow.launchDir}/results_ISMB_ess_and_delta_scan_long/${ess}_${delta}", mode: 'move'
	
	input:
	tuple val(ess), val(delta)
	
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
	config['Fitness']['Delta'] = ${delta}
	config['Signatures']['CNVs']['Events'][-1]['Prob'] = 0
	config['EvoParams']['PWGD'] = 0
	with open('config.json', 'w') as f:
		json.dump(config, f, indent=4)
	"
	~/.conda/envs/simcha/lib/dotnet/dotnet run --no-build --project ${simcha_path} -- -C config.json -D ${workflow.launchDir}/data/hg19 -e -R 5000 -O "." --light
        """
}

workflow {
	def ess = Channel.from(params.ess).take(50)
	def delta = Channel.from(params.delta).take(50)
	def product = ess.combine(delta)
	SimChA(product)
}
