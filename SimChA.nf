#!/usr/bin/env nextflow
nextflow.enable.dsl=2

simcha_path = workflow.launchDir + "/SimChA"

params.simcha_params_file = "default_params.json"

import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

def currentDateTime = LocalDateTime.now()
def formatter = DateTimeFormatter.ofPattern("yy_MM_dd_HH_mm_ss")
def timestamp = currentDateTime.format(formatter)

process SimChA {
    publishDir "${workflow.launchDir}/results/${timestamp}/${whole_chr_p}_${internal_p}_${telomere_p}", mode: 'move'

    input:
        val config
        tuple val(whole_chr_p), val(internal_p), val(telomere_p)

    output:
        path("*")

    script:
        def new_config = config        
        new_config.Signatures.each { sig ->
            // Iterate through the events
            sig.Events.each { event ->
                if(event.Type in ['ChromDeletion', 'ChromDuplication']) {
                    event.Prob = whole_chr_p
                }
                else if (event.Type in ['InternalDuplication', 'InternalDeletion', 'InternalInversion']) {
                    event.Prob = internal_p
                }
                else if (event.Type in ['TelomereDeletion', 'TelomereDuplication']) {
                    event.Prob = telomere_p
                }
            }
        }   
        def config_json = JsonOutput.toJson(new_config)
        """
        pwd
        echo '${config_json}' > config.json
        dotnet run --no-build --project ${simcha_path} -- -C config.json --data ${workflow.launchDir}/data -O "." -R 10
        """
}

workflow {
    params_file = file('default_params.json')
    def config = new JsonSlurper().parseText(params_file.text)
    
    def whole_chr_p = Channel.from(params.whole_chr_p)
    def internal_p = Channel.from(params.internal_p)
    def telomere_p = Channel.from(params.telomere_p)
    // Create a product of the above channels
    def product = whole_chr_p.combine(internal_p).combine(telomere_p)
    SimChA(config, product)
}
