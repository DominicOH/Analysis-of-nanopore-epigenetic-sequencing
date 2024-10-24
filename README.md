# Double and single stranded detection of 5-methylcytosine and 5-hydroxymethylcytosine with nanopore sequencing

These scripts and notebooks contain code used in the analysis of nanopore sequence data for the publication: Double and single stranded detection of 5-methylcytosine and 5-hydroxymethylcytosine with nanopore sequencing (_in press_). 

New nanopore sequence data in raw fast5 format as well as aligned bam format (aligned to mm39) has been made available on the NCBI SRA under the BioProject: PRJNA1144670. BedMethyl format data has also been uploaded to the NCBI GEO archive. 

## How to use this repository

This repository is intended to enable reproduction of figures and statistics used in the associated publication. Users are recommended to download the extracted modified base information from the NCBI GEO repository where they are archived as GSE279860. Users will then need to replace filepaths used in these scripts to match their filesystem. Scripts are named according to the statistics, figure, or supplementary figure that they produce or are related to.

Some reference and data files, excluding large genomic references for which links are provided, in `feature_references` and `data` folders. Where this is not possible due to size, instructions are provided within scripts. 

### If downloading raw fast5

Fast5 data requires conversion to pod5 format file using pod5 tools before base-calling and aligning (to the mm39/GRCm39 reference genome) with Dorado. Do so using the pod5-file-format package: https://github.com/nanoporetech/pod5-file-format. 

### If downloading bam format or aligning from basecalling

Pre-processing is required first by sorting according to chromosomal position (`samtools sort`) and indexing using `samtools index`.
 
### Extracting modified bases

Extracting modified bases was performed using `modkit pileup` and `modkit pileup-hemi` using the `--cpg` and `--mask` options, using the UCSC mm39 reference genome as a reference. 

## Comparison with public datasets

This repository refers to public datasets, including oxBS-seq and TAB-seq data produced by Ma _et al._ (2017). This was downloaded from the Beijing Genome Sequence Archive (GSA), where it is stored under the experimental accession numbers: CRX008031 (https://ngdc.cncb.ac.cn/gsa/browse/CRA000145/CRR008807) and CRX008030 (https://ngdc.cncb.ac.cn/gsa/browse/CRA000145/CRR008808). These datasets were processed using the following pipeline:

1. Raw data was downloaded in fastq.gz format.
2. Raw data was trimmed using Trim Galore! (see docs. https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
3. Trimmed fastq were aligned to a bisulphite converted reference genome (mm39) using `bismark` (see docs. https://www.bioinformatics.babraham.ac.uk/projects/bismark/). 
    -   Bisulphite converted reference genome prepared using `bismark_genome_preparation` 
4. Duplicate reads were removed using `deduplicate_bismark`
5. Modified bases were extracted in CpG positions using `bismark_methylation_extractor --merge_non_CpG --bedGraph --zero_based` 
6. Repeat sequences (repeatMasker) were removed using `bedtools intersect -v mm39_repeatMasker.bed`.   
