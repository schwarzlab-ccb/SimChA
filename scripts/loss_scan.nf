#!/usr/bin/env nextflow
nextflow.enable.dsl=2
simcha_path = workflow.launchDir + "/SimChA"
params.simcha_params_file = "/projects/ag-schwarzr/project-simcha/simcha/configs/losses_gains_separate.json"

import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

process SimChA {
	publishDir "${workflow.launchDir}/results_loss/${r_losses}_${alpha}", mode: 'move'
	
	input:
	val config
	tuple val(r_losses), val(alpha)
	
	output:
	path("*")
	
	script:
	def new_config = config
	new_config.Signatures.Losses.Prob = r_losses
	new_config.Fitness.TotalStrength = alpha
	new_config.Fitness.Stress = 1
        new_config.Fitness.TsgOg = 0
        new_config.Fitness.Essentiality = 0
	new_config.EvoParams.ThetaFitness = 20
	
	def config_json = JsonOutput.toJson(new_config)
	"""
	echo '${config_json}' > config.json
	dotnet run --no-build --project ${simcha_path} -- -C config.json -D ${workflow.launchDir}/data/hg19 -e -R 1000 -O "."
	"""
}

workflow {
	def params_file = file(params.simcha_params_file)
	def max_vals_per_param = 50
	def config = new JsonSlurper().parseText(params_file.text)
	println config
	def r_losses = Channel.from(params.r_losses).take(max_vals_per_param)
	def alpha = Channel.from(params.alpha).take(max_vals_per_param)
	def product = r_losses.combine(alpha)
	SimChA(config, product)
}
