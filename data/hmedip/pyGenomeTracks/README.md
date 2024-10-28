# Introduction
This section of the repository refers to how to reproduce the pyGenomeTracks image used in Figure 3a. 

# Data availability

Where possible, data has been provided in the `data/` directory; however, several files are not provided due to file size. To reproduce these, follow the steps below:
1.  Produce 5hmC bedgraph files using: `modkit pileup <in.bam> <out.bed> --bedgraph --cpg --ref <reference_genome.fa> --region chr7:45,746,000-45,751,000 --ignore m`
2. Produce depth .bg files using `samtools depth <in.bam> -r chr7:45,746,000-45,751,000 -d > <out.bg>`
3. Download `mm39.ncbiRefSeq.gtf.gz` from https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/genes/mm39.ncbiRefSeq.gtf.gz. 