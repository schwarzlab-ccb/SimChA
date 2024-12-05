#!/usr/bin/env nextflow
nextflow.enable.dsl=2
simcha_path = workflow.launchDir + "/SimChA"
params.simcha_params_file = "/projects/ag-schwarzr/project-simcha/simcha/configs/ess_scan.json"

import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

process SimChA {
	publishDir "${workflow.launchDir}/results_ess_scan/${ess}", mode: 'move'
	
	input:
	val config
	val ess
	
	output:
	path("*")
	
	script:
	a = 0.1
	delta = a/(1.0-ess)
	
	def new_config = config
	new_config.Fitness.TotalStrength = delta
	new_config.Fitness.Stress = 1.0-ess
        new_config.Fitness.TsgOg = 0
        new_config.Fitness.Essentiality = ess
	new_config.EvoParams.ThetaFitness = 20
	//new_config.Fitness.Haploinsufficiency = false
	def config_json = JsonOutput.toJson(new_config)
	"""
	echo '${config_json}' > config.json
	dotnet run --no-build --project ${simcha_path} -- -C config.json -D ${workflow.launchDir}/data/hg19 -e -R 4000 -O "."
	"""
}

workflow {
	def params_file = file(params.simcha_params_file)
	def max_vals_per_param = 170
	def config = new JsonSlurper().parseText(params_file.text)
	println config
	def ess = Channel.from(params.ess_1).take(max_vals_per_param)
	SimChA(config, ess)
}
