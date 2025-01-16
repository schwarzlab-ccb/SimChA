#!/usr/bin/env nextflow
nextflow.enable.dsl=2
simcha_path = workflow.launchDir + "/SimChA"
params.simcha_params_file = "/projects/ag-schwarzr/project-simcha/simcha/configs/tetraploid_rloss.json"

process SimChA {
	publishDir "${workflow.launchDir}/results_diploid_ess_scan/${ess}", mode: 'move'
	
	input:
	val ess
	
	output:
	path("*")
	
	script:
	"""
	cp ${params.simcha_params_file} config.json
	python -c "
	import json 
	with open('config.json', 'r') as f:
		config = json.load(f)
	config['Fitness']['Stress'] = 1.680
	config['Fitness']['TsgOg'] = 0
	config['Fitness']['Essentiality'] = ${ess}
	config['Fitness']['Haploinsufficiency'] = False
	config['Signatures']['CNVs']['Events'][0]['Prob'] = 4.5
	config['Signatures']['CNVs']['Events'][1]['Prob'] = 4.5
	config['Signatures']['CNVs']['Events'][2]['Prob'] = 4.5
	config['EvoParams']['TetraploidStart'] = False
	config['EvoParams']['EvolveInTime'] = False
	config['EventCountMean'] = 9
	with open('config.json', 'w') as f:
		json.dump(config, f, indent=4)
	"
	dotnet run --no-build --project ${simcha_path} -- -C config.json -D ${workflow.launchDir}/data/hg19 -e -R 10000 -O "." --light
        """
}

workflow {
	def max_vals_per_param = 1000
	def ess = Channel.from(params.ess).take(max_vals_per_param)
	SimChA(ess)
}
