#!/usr/bin/env nextflow


simcha_path = workflow.launchDir + "/SimChA"

params.simcha_params_file = "default_params.json"

params.assembly = ["hg19", "hg38"]

import groovy.json.JsonSlurper
import groovy.json.JsonOutput

process SimChA {
    input:
        val config
        val assembly

    output:
        path("config.json")

    script:
        def new_config = config
        new_config.Assembly = assembly
        def config_json = JsonOutput.toJson(new_config)
        """
        pwd
        echo '${config_json}' > config.json
        dotnet run --no-build --project ${simcha_path} -- -C config.json --data ${workflow.launchDir}/data
        """
}

workflow {
    params_file = file('default_params.json')
    def config = new JsonSlurper().parseText(params_file.text)
    SimChA(config, Channel.from(params.assembly))
}
