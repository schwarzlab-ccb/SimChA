#!/usr/bin/env nextflow

simcha_path = workflow.launchDir + "/SimChA"

params.simcha_params_file = "default_params.json"

params.assembly = ["hg19", "hg38"]

import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import java.time.LocalDateTime
import java.time.format.DateTimeFormatter

def currentDateTime = LocalDateTime.now()
def formatter = DateTimeFormatter.ofPattern("yy_MM_dd_HH_mm_ss")
def timestamp = currentDateTime.format(formatter)

process Build {
    script:
        """
        dotnet build ${simcha_path}
        """
    output:
        stdout
}

process SimChA {
    publishDir "${workflow.launchDir}/results/${timestamp}/${assembly}", mode: 'copy'

    input:
        val build_out
        val config
        val assembly

    output:
        path("*")

    script:
        def new_config = config
        new_config.Assembly = assembly
        def config_json = JsonOutput.toJson(new_config)
        """
        pwd
        echo '${config_json}' > config.json
        dotnet run --no-build --project ${simcha_path} -- -C config.json --data ${workflow.launchDir}/data -O "."
        """
}

workflow {
    def build_out = Build()
    params_file = file('default_params.json')
    def config = new JsonSlurper().parseText(params_file.text)
    SimChA(build_out, config, Channel.from(params.assembly))
}
