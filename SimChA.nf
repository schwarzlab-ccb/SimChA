#!/usr/bin/env nextflow
nextflow.enable.dsl=2

simcha_path = workflow.launchDir + "/SimChA"
params.simcha_params_file =  "default_params.json"

import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

def currentDateTime = LocalDateTime.now()
def formatter = DateTimeFormatter.ofPattern("yy_MM_dd_HH_mm_ss")
def timestamp = currentDateTime.format(formatter)

process SimChA {
    publishDir "${workflow.launchDir}/results/${timestamp}/${chrDup}_${chrDel}_${intDup}_${intDel}_${tailDel}_${bfb}_${wgd}", mode: 'move'

    input:
        val params_text
        tuple val(chrDup), val(chrDel), val(intDup), val(intDel), val(tailDel), val(bfb), val(wgd)

    output:
        tuple path("*")

    script:
        def new_config = new JsonSlurper().parseText(params_text)     
        new_config.Signatures.each { key, val ->
            // Iterate through the events
            val.Events.each { event ->
                if(event.Type in ['ChromDeletion']) {
                    event.Prob = chrDel
                }
				else if (event.Type in ['ChromDuplication']) {
					event.Prob = chrDup
				} 
				else if (event.Type in ['InternalDeletion']) {
					event.Prob = intDel
				}
				else if (event.Type in ['InternalDuplication']) {
					event.Prob = intDup
				}
                else if (event.Type in ['TailDeletion']) {
                    event.Prob = tailDel
                }
				else if (event.Type in ['BreakageFusionBridge']) {
                    event.Prob = bfb
                }
                else if (event.Type in ['WholeGenomeDoubling']) {
                    event.Prob = wgd
                }
                else {
                    event.Prob = 0
                }
            }
        }   
        def config_json = JsonOutput.toJson(new_config)
        """
        echo '${config_json}' > config.json
        dotnet run --no-build --project ${simcha_path} -- -C config.json --data ${workflow.launchDir}/data/hg19 -R ${params.reps} -O "."
        """
}

workflow {
    def params_file = file('default_params.json')
    def max_params = 2
    def chrDup = Channel.from(params.chrDup).take(max_params)
	def chrDel = Channel.from(params.chrDel).take(max_params)
	def intDup = Channel.from(params.intDup).take(max_params)
	def intDel = Channel.from(params.intDel).take(max_params)
	def tailDel = Channel.from(params.tailDel).take(max_params)
	def bfb = Channel.from(params.bfb).take(max_params)
	def wgd = Channel.from(params.wgd).take(max_params)
    def product = chrDup
		.combine(chrDel)
		.combine(intDup)
		.combine(intDel)
		.combine(tailDel)
		.combine(bfb)
		.combine(wgd)
    SimChA(params_file.text, product)
}
