#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process Test {
    input:
        val my_val
        
    output:
        stdout 

    script:
        """
        echo ${my_val}
        """
}

workflow {
    // Assuming channels A, B, C are defined
    channelA = Channel.from(1, 2, 3)
    channelB = Channel.from('a', 'b', 'c')
    channelC = Channel.from('X', 'Y', 'Z')

    // Cross channelA and channelB
    crossAB = channelA.combine(channelB)

    // Then cross the result with channelC
    result = crossAB.combine(channelC)

    // Consume the result
    result.view()

    Test(result)
}