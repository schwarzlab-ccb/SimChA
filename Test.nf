#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process Test {        
    output:
        file("*")

    script:
        """       
        echo $workflow.launchDir > test.echo
        echo $workflow.launchDir > test.1
        echo $workflow.launchDir > test.2
        echo $workflow.launchDir > test.3
        echo $workflow.launchDir > test.4
        echo $workflow.launchDir > test.5
        echo $workflow.launchDir > test.6
        dotnet run --project ${workflow.launchDir}/TestApp  -- -C config.json --data
        """
}

workflow {
    Test().view()
}