#!/usr/bin/env nextflow
nextflow.enable.dsl=2
simcha_path = workflow.launchDir + "/SimChA"
params.simcha_params_file = "/projects/ag-schwarzr/project-simcha/simcha/configs/wgd_and_loss.json"

process SimChA {
	publishDir "${workflow.launchDir}/scripts/results_roland_exp_tailless/${r_loss}_${alpha}_${r_wgd}_${mut_rate}", mode: 'move'
	
	input:
	tuple val(r_loss), val(alpha), val(r_wgd), val(mut_rate)
	
	output:
	path("*")
	
	script:
	ess = 1.0 - alpha
	
	"""
	cp ${params.simcha_params_file} config.json
	python -c "
	import json
	with open('config.json', 'r') as f:
		config = json.load(f)
	config['Fitness']['Stress'] = ${alpha}
	config['Fitness']['TsgOg'] = 0
	config['Fitness']['Essentiality'] = ${ess}
	config['Fitness']['TotalStrength'] = 1
	config['Fitness']['Haploinsufficiency'] = True
	config['Signatures']['CNVs']['Events'][0]['Prob'] = ${r_loss}
	config['Signatures']['CNVs']['Events'][1]['Prob'] = ${r_loss}
	config['Signatures']['CNVs']['Events'][-1]['Prob'] = ${r_wgd}
	config['EvoParams']['MutationRate'] = ${mut_rate}
	with open('config.json', 'w') as f:
		json.dump(config, f, indent=4)
	"
	dotnet run --no-build --project ${simcha_path} -- -C config.json -D ${workflow.launchDir}/data/hg19 -e -R 2000 -O "." --light
	"""
}

workflow {
	def max_vals_per_param = 20
	def r_loss = Channel.from(params.r_loss).take(max_vals_per_param)
	def alpha = Channel.from(params.alpha).take(max_vals_per_param)
	def delta = Channel.from(params.delta).take(max_vals_per_param)
	def product = r_loss.combine(alpha).combine(r_wgd).combine(mut_rate)
	SimChA(product)
}
