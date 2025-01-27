# Introduction
This section of the repository refers to how to reproduce the pyGenomeTracks image used in Figure 3a. 

# Data availability

Where possible, data has been provided in `fig3a_source.tar.gz` (see Source Data on the front page). The file `wgs_ref_h` is not provided in this archive, but can be created as follows: 
1. Download the input bam file from Sequence Read Archive accession SRR30150148 (11.4 Gb). 
2. Sort and index the bam using `samtools sort <input_bam> && samtools index <sorted_input_bam>` 
2. Produce 5hmC bedgraph using: `modkit pileup <sorted.bam> <out.bed> --bedgraph <output directory> --cpg --ref <reference_genome.fa> --ignore m`
3. (Optional) limit the size of the output to the selected region using: `--region chr7:45,746,000-45,751,000`

Additionally, the genome reference file: `mm39.ncbiRefSeq.gtf.gz` requires downloading from: https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/mm39.ncbiRefSeq.gtf.gz. 