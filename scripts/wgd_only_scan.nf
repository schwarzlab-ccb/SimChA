#!/usr/bin/env nextflow
nextflow.enable.dsl=2
simcha_path = workflow.launchDir + "/SimChA"
params.simcha_params_file =  "/projects/ag-schwarzr/project-simcha/simcha/configs/wgd_only.json"

process SimChA {
	publishDir "${workflow.launchDir}/results_wgd_and_alpha/${mut_rate}_${alpha}", mode: 'move'
	
	input:
	tuple val(mut_rate), val(alpha)
	
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
	config['EvoParams']['MutationRate'] = ${mut_rate}
	config['EvoParams']['EvolveInTime'] = True
	config['EvoParams']['MaxTime'] = 100
	config['EvoParams']['WGDAccelerationFactor'] = 4
	with open('config.json', 'w') as f:
		json.dump(config, f, indent=4)
        "
        dotnet run --no-build --project ${simcha_path} -- -C config.json -D ${workflow.launchDir}/data/hg19 -e -R 2000 -O "." --light
        """
}

workflow {
	def max_vals_per_param = 50
	def mut_rate = Channel.from(params.mut_rate).take(max_vals_per_param)
	def alpha = Channel.from(params.alpha).take(max_vals_per_param)
	def product = mut_rate.combine(alpha)
	SimChA(product)
}
