#!/usr/bin/env nextflow
nextflow.enable.dsl=2
simcha_path = workflow.launchDir + "/SimChA"
params.simcha_params_file = "/projects/ag-schwarzr/project-simcha/simcha/configs/tetraploid_arm_and_chrom.json"

import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

process SimChA {
	publishDir "${workflow.launchDir}/results_tetraploid_arm_and_chrom/${r_loss}_${alpha}", mode: 'move'
	
	input:
	val config
	tuple val(r_loss), val(alpha)
	
	output:
	path("*")
	
	script:
	def new_config = config
	// new_config.Signatures.Losses.Prob = r_loss
	// new_config.Signatures.Gains.Prob = 1
	new_config.Signatures.CNVs.Events.each { event ->
		if (event.Type in ['ChromDeletion', 'ArmDeletion']) {
			event.Prob = r_loss
		}
		if (event.Type in ['ChromDuplication', 'ArmDuplication']) {
			event.Prob = 1
		}
	}
	new_config.Fitness.TotalStrength = alpha
	new_config.Fitness.Stress = 1
        new_config.Fitness.TsgOg = 0
        new_config.Fitness.Essentiality = 0
	new_config.EvoParams.ThetaFitness = 20
	
	def config_json = JsonOutput.toJson(new_config)
	"""
	echo '${config_json}' > config.json
	dotnet run --no-build --project ${simcha_path} -- -C config.json -D ${workflow.launchDir}/data/hg19 -e -R 15000 -O "."
	"""
}

workflow {
	def params_file = file(params.simcha_params_file)
	def max_vals_per_param = 35
	def config = new JsonSlurper().parseText(params_file.text)
	println config
	def r_loss = Channel.from(params.r_loss_long).take(max_vals_per_param)
	def alpha = Channel.from(params.alpha_long).take(max_vals_per_param)
	def product = r_loss.combine(alpha)
	SimChA(config, product)
}
