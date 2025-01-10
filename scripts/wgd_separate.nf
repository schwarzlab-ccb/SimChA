#!/usr/bin/env nextflow
nextflow.enable.dsl=2
simcha_path = workflow.launchDir + "/SimChA"
params.simcha_params_file = "/projects/ag-schwarzr/project-simcha/simcha/configs/wgd_separated_events.json"

import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

process SimChA {
	publishDir "${workflow.launchDir}/results_wgd_separated_events_mut_1/${ess}", mode: 'move'
	
	input:
	val config
	val ess
	
	output:
	path("*")
	
	script:
	println ess
	
	def new_config = config
	new_config.Fitness.TotalStrength = 1
	new_config.Fitness.TsgOg = 1.0-ess
        new_config.Fitness.Stress = 0
        new_config.Fitness.Essentiality = ess
	new_config.EvoParams.ThetaFitness = 20
	new_config.EvoParams.MutationRate = 1
	new_config.Fitness.Haploinsufficiency = true
	def config_json = JsonOutput.toJson(new_config)
	
	"""
	echo '${config_json}' > config.json
	dotnet run --no-build --project ${simcha_path} -- -C config.json -D ${workflow.launchDir}/data/hg19 -e -R 10000 -O "." --light
	"""
}

workflow {
	def params_file = file(params.simcha_params_file)
	def max_vals_per_param = 1
	def config = new JsonSlurper().parseText(params_file.text)
	println JsonOutput.prettyPrint(JsonOutput.toJson(config))
	def ess = Channel.from(params.ess_1).take(max_vals_per_param)
	println "params.ess: ${params.ess}"
	SimChA(config, ess)
}
