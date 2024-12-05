#!/usr/bin/env nextflow
nextflow.enable.dsl=2
simcha_path = workflow.launchDir + "/SimChA"
params.simcha_params_file = "/projects/ag-schwarzr/project-simcha/simcha/configs/wgd_and_loss.json"

import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

process SimChA {
	publishDir "${workflow.launchDir}/results_wgd_and_loss/${r_loss}_${alpha}_${ess}_${mut_rate}", mode: 'move'
	
	input:
	val config
	tuple val(r_loss), val(alpha), val(ess), val(mut_rate)
	
	output:
	path("*")
	
	script:
	def new_config = config
	new_config.Signatures.CNVs.Events.each { event ->
		if (event.Type in ['WholeGenomeDoubling']) {
			event.Prob = 0.5
		}
		if (event.Type in ['ChromDeletion']) {
			event.Prob = r_loss
		}
	}
	new_config.Fitness.TotalStrength = 1
	new_config.Fitness.Stress = alpha
        new_config.Fitness.TsgOg = 0
        new_config.Fitness.Essentiality = ess
	new_config.EvoParams.ThetaFitness = 20
	new_config.EvoParams.MutationRate = mut_rate
	
	def config_json = JsonOutput.toJson(new_config)
	"""
	echo '${config_json}' > config.json
	dotnet run --no-build --project ${simcha_path} -- -C config.json -D ${workflow.launchDir}/data/hg19 -e -R 500 -O "."
	"""
}

workflow {
	def params_file = file(params.simcha_params_file)
	def max_vals_per_param = 10
	def config = new JsonSlurper().parseText(params_file.text)
	println config
	def r_loss = Channel.from(params.r_loss).take(max_vals_per_param)
	def alpha = Channel.from(params.alpha).take(max_vals_per_param)
	def ess = Channel.from(params.ess).take(max_vals_per_param)
	def mut_rate = Channel.from(params.mut_rate).take(max_vals_per_param)
	def product = r_loss.combine(alpha).combine(ess).combine(mut_rate)
	SimChA(config, product)
}
