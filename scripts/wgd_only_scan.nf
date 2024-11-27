#!/usr/bin/env nextflow
nextflow.enable.dsl=2
simcha_path = workflow.launchDir + "/SimChA"
params.simcha_params_file =  "/projects/ag-schwarzr/project-simcha/simcha/configs/wgd_only.json"

import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

process SimChA {
	publishDir "${workflow.launchDir}/results_wgd/${mut_rate}_${alpha}", mode: 'move'
	
	input:
	val config
	tuple val(mut_rate), val(alpha)
	
	output:
	path("*")
	
	script:
	def new_config = config
	new_config.Fitness.TotalStrength = alpha
	new_config.Fitness.Stress = 1
	new_config.Fitness.TsgOg = 0
	new_config.Fitness.Essentiality = 0
	new_config.EvoParams.MutationRate = mut_rate
	
	def config_json = JsonOutput.toJson(new_config)
	"""
	echo '${config_json}' > config.json
	dotnet run --no-build --project ${simcha_path} -- -C config.json -D ${workflow.launchDir}/data/hg19 -e -R 500 -O "."
	"""
}

workflow {
	def params_file = file(params.simcha_params_file)
	def max_vals_per_param = 50
	def config = new JsonSlurper().parseText(params_file.text)
	println config
	def mut_rate = Channel.from(params.mut_rate).take(max_vals_per_param)
	def alpha = Channel.from(params.alpha).take(max_vals_per_param)
	def product = mut_rate.combine(alpha)
	SimChA(config, product)
}
