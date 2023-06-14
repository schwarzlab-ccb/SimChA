#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process Test {
    input:
        val my_val
        
    output:
        stdout 

    script:
        """
        dotnet -h
        """
}

workflow {
    Test("test").view()
}