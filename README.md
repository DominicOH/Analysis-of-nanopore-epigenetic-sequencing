# Double and single stranded detection of 5-methylcytosine and 5-hydroxymethylcytosine with nanopore sequencing

These scripts and notebooks contain code used in the analysis of nanopore sequence data for the publication: Double and single stranded detection of 5-methylcytosine and 5-hydroxymethylcytosine with nanopore sequencing (_in press_). 

New nanopore sequence data in raw fast5 format as well as aligned bam format (aligned to mm39) has been made available on the NCBI SRA under the BioProject: PRJNA1144670. BedMethyl format data has also been uploaded to the NCBI GEO archive. 

## How to use this repository

This repository is intended to aid reproduction of figures and statistics used in the associated publication. Users are recommended to download the extracted modified base information from the NCBI GEO repository where they are archived as GSE279860. Users will then need to replace filepaths used in these scripts to match their filesystem. Scripts are named according to the statistics, figure, or supplementary figure that they produce or are related to.

With the exception of large genomic references, for which links are provided, references and data files are provided in `feature_references` and `data` folders. Where this is not possible due to size, instructions are provided within scripts. 

## Where to find data

Follow the README.md file in the `data/` directory for links and information to access all datasets. 

## Comparison with public datasets

This repository refers to public datasets, including oxBS-seq and TAB-seq data produced by Ma _et al._ (2017). This was downloaded from the Beijing Genome Sequence Archive (GSA), where it is stored under the experimental accession numbers: CRX008031 (https://ngdc.cncb.ac.cn/gsa/browse/CRA000145/CRR008807) and CRX008030 (https://ngdc.cncb.ac.cn/gsa/browse/CRA000145/CRR008808). These datasets were processed using the following pipeline:

1. Raw data was downloaded in fastq.gz format.
2. Raw data was trimmed using Trim Galore! (see docs. https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
3. Trimmed fastq were aligned to a bisulphite converted reference genome (mm39) using `bismark` (see docs. https://www.bioinformatics.babraham.ac.uk/projects/bismark/). 
    -   Bisulphite converted reference genome prepared using `bismark_genome_preparation` 
4. Duplicate reads were removed using `deduplicate_bismark`
5. Modified bases were extracted in CpG positions using `bismark_methylation_extractor --merge_non_CpG --bedGraph --zero_based` 
6. Repeat sequences (repeatMasker) were removed using `bedtools intersect -v mm39_repeatMasker.bed`.   

Additionally, we make use of hMeDIP-seq data procured from NCBI's Gene Expression Omnibus under accession: GSE25398. These datasets were processed using a pipeline comprised of Trim Galore!, Picard, and MACS2, which is available here: https://github.com/DominicOH/ChIP2MACS2. 
