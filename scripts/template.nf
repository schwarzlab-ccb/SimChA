#!/usr/bin/env nextflow
nextflow.enable.dsl=2

simcha_path = workflow.launchDir + "/SimChA"
params.simcha_params_file =  workflow.launchDir + "/configs/losses_gains_separate.json"
w_wgd = 0.35
w_total_loss = 2.73 + 2.47 + 5.93 + 29.6
w_total_gain = 2.12 + 1.9 + 6.01 + 49.2

import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

process SimChA {
	publishDir "${workflow.launchDir}/results/${losses_p}_${ess}", mode: 'move'"
	
	input:
	val config
	tuple val(losses_p), val(ess)
	
	output:
	path("*")
	
	script:
	def w_total = w_total_loss * losses_p + w_total_gain
	def new_config = config
	new_config.Signatures.Losses.Prob = losses_p
	new_config.Signatures.WGD.Prob = w_wgd * w_total / (1.0 - w_wgd)
	new_config.Fitness.Essentiality = ess
	
	def config_json = JsonOutput.toJson(newconfig)
	"""
	echo '${config_json}' > config.json
	dotnet run --no-build --project ${simcha_path} -- -C config.json -D ${workflow.launchDir}/data/hg19 -e -R 2000 -O "."
	"""
}

workflow {
	def params_file = file(params.simcha_params_file)
	def max_vals_per_param = 1
	def config = new JsonSlurper().parseText(params_file.text)
	def losses_p = Channel.from(params.losses_p).take(max_vals_per_param)
	def ess = Channel.from(params.ess).take(max_vals_per_param)
	def product = losses_p.combine(ess)
	SimChA(config, product)
}
