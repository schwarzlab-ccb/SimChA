#!/usr/bin/env nextflow
nextflow.enable.dsl=2
simcha_path = workflow.launchDir + "/SimChA"
params.simcha_params_file = "/projects/ag-schwarzr/project-simcha/simcha/configs/tetraploid_rloss.json"

process SimChA {
	publishDir "${workflow.launchDir}/results_tetraploid_rloss/${r_loss}_${alpha}", mode: 'move'
	
	input:
	tuple val(r_loss), val(alpha)
	
	output:
	path("*")
	
	script:

	"""
	cp ${params.simcha_params_file} config.json
	python -c "
	import json
	with open('config.json', 'r') as f:
		config = json.load(f)
	config['Fitness']['Stress'] = ${alpha}
	config['Fitness']['TsgOg'] = 0
	config['Fitness']['Essentiality'] = 0
	config['Fitness']['Haploinsufficiency'] = False
	config['Signatures']['CNVs']['Events'][0]['Prob'] = ${r_loss}
	config['Signatures']['CNVs']['Events'][1]['Prob'] = ${r_loss}
	config['Signatures']['CNVs']['Events'][2]['Prob'] = ${r_loss}
	config['EvoParams']['TetraploidStart'] = True
	config['EvoParams']['EvolveInTime'] = False
	config['EventCountMean'] = 35
	with open('config.json', 'w') as f:
		json.dump(config, f, indent=4)
	"
	dotnet run --no-build --project ${simcha_path} -- -C config.json -D ${workflow.launchDir}/data/hg19 -e -R 10000 -O "." --light
        """
}

workflow {
	def max_vals_per_param = 50
	def r_loss = Channel.from(params.r_loss).take(max_vals_per_param)
	def alpha = Channel.from(params.alpha).take(max_vals_per_param)
	def product = r_loss.combine(alpha)
	SimChA(product)
}
