# SimChA: Simulator of Chromosomal Aberrations

SimChA is a research project at SchwarzLab.

## Requirements

### TL;DR; 

The program can be run on a platform of your choice in the provided *Conda* environment. Inside of the repo run 
```
conda env create --file simcha.yml
conda activate simcha
```

### Tested platforms

The program has been tested on:
* Windows 10 - PowerShell
* Windows 10 - WSL2 Ubuntu
* Ubuntu 20
* MacOS X 10 

### Gene Data

The gene data (e.g. PCAWG, TSGs and OGs) is stored using the Git LFS system, so in order for SimChA to work, you will need to download Git LFS (https://git-lfs.com/) and once you have install Git LFS, inside the repo, run:
```
git lfs pull
```

### Somatic Nucleotide Variants (SNVs)

SimChA is capable of handling the Jukes-Cantor nucleotide substitution model, but to attach proposed SNV events from SimChA to the reference genome, the reference genomes have to be downloaded. To do this, we have provided the `DownloadRefData.sh` script, which can be run as follows:

```
chmod +x DownloadRefData.sh && ./DownloadRefData.sh
```

Note that if you don't want to download both reference genomes (hg19 and hg38), simply comment out or remove the relevant section of the script.


### Simulation

The simulation code is written in **C# 8**. **.NET 6.1** is requred. We recommend installation using *Conda*:
```
conda install -c conda-forge dotnet
``` 

### Data analysis

The analysis code is written in **Python 3.8**. The following packages are required (either from *Conda* or *Pip*):
```
conda install -c conda-forge -c bioconda biopython matplotlib numpy pandas seaborn pillow lifelines
```

## Execution

The default execution is:
```
git clone git@bitbucket.org:schwarzlab/simcha.git
cd simcha
dotnet run --project SimChA
./plots.sh
```

The results will be written to the folder `./out`

### Options

Use `dotnet run -- [options]` to specify any of the following:

```
  -O, --output    (Default: ./out) The path to the output files.
  -C, --config    (Default: ./default_params.json) A json file with configuration of the experiment.
  -T --tree path     A path to a newick with a phylogenetic tree. The file must contain columns "ID", "ParentID", and "Distance". For MCMC, also "Fitness" must be included.
  -R --repeats n   Will run the simulation and generate n independent samples.
  -M --mcmc      Will run the MCMC algorithm to match the clode fitness to the fitness given in the tree file or 1 if no tree file is given.
  -P --profile path   A path to a file with copy number profiles.
  -D --data (Default: ./data/hg19) The path for the tsv-files with the chromosome list, and the essential-, tsg- and og-genes with scores. The files should be named: chromosomes.tsv, essential.tsv, tsg.tsv, og.tsv and contained in a folder with the same name as the assembly used. If you want to use either of the flags "--variants" or "--fasta", you must also have the corresponding reference genome, named genome.fa in the same directory.
  -s  Calculate consistent copy numbers segmentation in addition to a normal CN file. The output file, consistent_CNs.tsv, will have NA if the original sample did not have data in a given region. Note! This file is going to be very large if you have many samples.
  -V, --variants     Generate a Variant Call Format-like (VCF) file at the end of the run.
  -F, --fasta       Generate the FASTA file of the produced karyotype (Note: this currently only works for single-sample simulations).
```

### Parameters

Default parameters are found in the file: `default_params.json`. The exact parameters are dependent on the execution mode. The parameters are:

#### Common (R, T, P)
* `Fitness: { Stress: int, TsgOg: int, Essentiality: int }` - The fitness function to use. The fitness is calculated as the product of the three scores. The scores are weighted by the parameters.

#### Simulation mode only (R, T)
* `Seed: int` - The seed for the random number generator.
* `Sex: ["Female", "Male", "None"]` - Whether to simulate exclusively female chromosome, male chromosome, or randomly select from the two.
* `AutosomesOnly: bool` - Choose whether to simulate sex chromosomes or not.

#### Repeats mode only (R)
* `EventCount: int` - The mean number of events per sample.
* `EventDist: ["Uniform", "Exponential", "Normal"]` - The distribution of the number of events per sample. Uniform means all with be exactly the `EventCount`.

#### MCMC mode (M)
TODO

## Input Data

Data should be provided as a path to a folder using the `--data` parameter. The folder needs to contain description of the chromosomes and TSG, OG, and Essentiality scores for the individual genes. We provide the `GRCh38.p14` and `GRCh38.p13` in the folders `./data/hg19` and `.data/hg38` respectivelly. 

### Chromosomes

The chromosome contains the lengths and sex-specific information of individual chromosomes. By default there are two columns, one with the name of a chromosome, one with a number of basises.

Unless the third column is provided, the penultimate chromosome is considered to be the X chromosome and the ultimate to be Y. All the other are considers to be autosomes. If you want to specify the information directly, the third column can contain a string with on of the following `Female`, `Male`, `Both`.

Example file:

```
chr1	248956422	Both
...
chrY	59373566	Male
```

### TSG/OG/Eseentiality Score
We use three data files, providing the TSG/OG/Essentiality score. Each file is a tab-separated file with the following columns:

* `Chromosome: string` - The name of the chromosome.
* `Start: long` - The start position of the gene (inclusive).
* `End: long` - The end position of the gene (inclusive).
* `Gene: string` - The name of the gene.
* `Score: float` - The score of the gene.

Example file:
```
chr3	178865902	178957881	PIK3CA	1
chr7	140419127	140624564	BRAF	0.991919559
chr12	25357723	25403870	KRAS	0.990684897
...
```

For each gene, the:

* TSG score is the probability that the gene is a Tumor Suppresor Gene (TSG) 
* OG score is the probability that the gene is an Oncogene (OG) 
* Essentiality score is the log2 fold change in the reproducibility of a cell after a knock-out of the gene (both alleles).

## Signatures

SimChA simulates events based on mutational signatures. Each signature is a set of events and their associated parameters, for example consider the following excerpt from a configuration file:

```
"Signatures" : {
  "WChE_example" : {
    "Prob": 1,
    "Events": [
      {
        "Type": "ChromDeletion",
        "Prob": 1
      },
      {
        "Type": "ChromDuplication",
        "Prob": 2
      }        
    ]
  },
  "Internal_example" : {
    "Prob": 5,
    "Events": [
      {
        "Type": "InternalDeletion",
        "Prob": 1
        "Size": 1000000
      }      
    ]
  }
}
```
These are two signatures, one for Whole Chromosome Events (WChE) and one for Internal events. The likelihood of a signature being selected is 1 : 5, meaning the Internal is 5 times as likely. 

If the WChE is selected, the probability of a deletion compared to duplication is 1 : 2. SimChA works with *contigs*, i.e. contiguous sequences of bases. A deletion or duplication will be a deletion or duplication of a contig. A contig may be comprised of parts of different chromosomes, e.g. after a translocation, however at the start of simulation, the set of contigs is the same as the set of chromosomes.

If the Internal Deletion is selected, the mean size of a deleted segment will be 1MB, distributed exponentially.

Each Event has an associated type and a probability, which is relative to all the other events in the signature. The following events are available:

* `ChromDeletion`
* `ChromDuplication`
* `TailDeletion`
* `InternalDeletion`
* `InternalDuplication`
* `InternalInversion`
* `InvertedDuplication`
* `Translocation`
* `BreakageFusionBridge`
* `WholeGenomeDoubling`
* `Chromothripsis`
* `Chromoplexy`
* `TIChain`
* `TICycle`
* `TIBridge`
* `Pyrgo`
* `Rigma`
* `SNV`

Events have parameters from the following:
* `Type: string` - The type of the event, one from the list above.
* `Prob: double` - The probability of the event being selected.
* `Size: long` - The mean size of the event in number of bases, exponentially distributed. Only applicable to internal events.
* `Frag: double` - Some complex events cause fragmentation, this is the mean number of fragments.

### Chromosome Deletion (Prob)

A single contig is selected at random and removed.

### Chromosome Duplication (Prob)

A single contig is selected at random and duplicated.

### Whole Genome Doubling (Prob)

All the existing contigs are duplicated.

### Tail Deletion (Prob, Size)

A tail of a length given by the `Size` parameter is removed from an end of a contig. The end is selected by a coin flip.

### Internal Deletion (Prob, Size)

A single contig is selected, from which a segment distributed by along the `Size` parameter is removed. The position of a segments is guaranteed to be internal and uniformly distributed.


### Internal Duplication (Prob, Size)
A single contig is selected, from which a segment distributed by along the `Size` parameter is duplicated. The position of a segments is guaranteed to be internal and uniformly distributed. This segment is pasted directly after its original position.

![Translocation Cartoon](./doc/img/InternalDuplication.png)

### Internal Inversion (Prob, Size)
A single contig is selected, from which a segment distributed by along the `Size` parameter is inverted. The position of a segments is guaranteed to be internal and uniformly distributed.

### Inverted Duplication (Prob, Size)
Like a duplication, but the segment is inverted before being pasted.

### Translocation  (Prob, Size)
Two contigs are selected, from which a segment distributed along the `Size` parameter is swapped. The position is selected individualy for each contig. A coin flip decides if one of the segments is inverted before being pasted.

![Translocation Cartoon](./doc/img/translocation.png)

### BreakageFusionBridge (Prob, Size)

A contig and its tail is selected (see above). The tail is then removed, the rest is copied, the copy is reversed and the two copies are connected on the breakage location.

![Translocation Cartoon](./doc/img/BFB_cartoon.png)

### Chromothripsis (Prob, Size)

1. A contig is broken into a number of fragments, such that the fragment size is distributed exponentially with a mean of `Size`.
2. Have, `f` the number of fragments from step 1). A `k` framgents are then randomly selected such that `0 < k <= f`.
3. The `k` fragments are then reassembled in a random order.  

TODO: image

### Chromoplexy (Prob, Size)

1. A number `c` of contigs is selected, following the probability distribution listed below.
2. Each contig is broken into a number of fragments, such that the fragment size is distributed exponentially with a mean of `Size`.
3. These contigs are then reassembled in a random order into a single contig.
4. This contig is broken down a uniformly distributed positions back into `c` contigs.	

Probabilities are sourced from [Ashby et al., 2019](https://ashpublications.org/blood/article/134/Supplement_1/3767/424006/Chromoplexy-and-Chromothripsis-Are-Important):

* 2 contigs: 0.23%
* 3 contigs: 0.46%
* 4 contigs: 0.16%
* 5 contigs: 0.10%
* 6 contigs: 0.05%

TODO: check

### Pyrgo 
![Pyrgo](./doc/img/pyrgo_128.png)
TODO

### Rigma
![Rigma](./doc/img/rigma_128.png)
TODO

### Template Insertion Chain
TODO

### Template Insertion Cycle
TODO

### Template Insertion Bridge
TODO


## Output

The text files are primarily used as source for plots shown below.

#### `copynumbers.tsv`

CN output in the format examplifed as:


| sample_id | chrom | start | end   | cn_a | cn_b | n_snvs |
|-----------|-------|-------|-------|------|------|--------|
|0          |	chr1  |	24721 |	98434 |	1    |   	4 |      2 |

#### `clone_tree.new`

An evolutionary tree with mutation distances and population sizes between the individual subclones. The graph is written in the [newick format](https://en.wikipedia.org/wiki/Newick_format).


#### `sim_params.json` 

Stores configuration parameters used for this simulation, including the random seed. If this file is provided on input, the exact same simulation will be executed.

#### `clones.tsv`  

Information about the individual subclones at the end of the simulation. This file has the chromosome ranges. A range is in the format:

`ChromID*[start:end)`

where `*` is either `+` for 5' to 3' or `-` for 3' to 5'.

The start position is inclusive, the end exclusive.


#### `events.tsv`

Information about the individual events in the simulation. The columns are:

| Id | 	EventType	| EventString	| dFitness |	tFitness |	Mutations |
|----|------------|-------------|----------|-----------|------------|
|  0 | ChromDuplication | contig:23 | 11.53782375 | 11.53782375 | 1 |

#### 'vcf.tsv'

Generates a simple VCF file of all the observable (i.e. final and nucleotide-changing) SNVs that occurred to the samples during simulation. The columns are:

`##fileformat=VCFv4.3`

`##source=SimChAV1.0`

`##reference=verily_hg19_genome.fa`

| `#SAMPLEID` | 	CHROM	| POS	| ID |	REF |	ALT |
|----|------------|-------------|----------|-----------|------------|
|  0 | chr1 | 8114090 | . | a | C |

#### 'genome.fa'

*Note that this only currently works with single-sample simulations* 

Generates the final simulated karyotype sequence from the relevant human genome reference. It will also include 


### Plots (created by `./plots.sh`)

In order to run `./plots.sh` you first have to install the consistent segmentation routine with by running `python setup.py build_ext --inplace` in the `simcha/scripts/segmentation_cython` directory (with [Cython](https://anaconda.org/anaconda/cython) installed).

#### `copy_numbers_heatmap.png` 
The copy number tracks

## Contact
Email questions, feature requests and bug reports to Adam Streck, adam.streck@mdc-berlin.de.

## License
SimChA is available under the MIT License.
