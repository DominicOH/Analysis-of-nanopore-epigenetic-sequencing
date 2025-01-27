# Data
Several datasets are referred to throughout the analysis.

## Tabulated modification data
We recommend downloading these processed data files (Series GSE279860):

- CBM2_1: [GSM8582102](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8582102)
- CBM2_2: [GSM8582103](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8582103)
- CBM3_1: [GSM8582104](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8582104)
- CBM3_2: [GSM8582105](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8582105)

## Alignment files
Most analyses can be completed using only the Gene Expression Omnibus bedMethyl files; however, if necessary, these BAM files can be downloaded as well. 

### DNA Methylation Standards
These alignment files were produced by sequencing the Human Methylated & Non-Methylated (WGA) DNA Set (Zymo Research, D5013).

- ZU1: [SRX26527543](https://www.ncbi.nlm.nih.gov/sra/SRX26527543[accn])
- ZU2: [SRX26527544](https://www.ncbi.nlm.nih.gov/sra/SRX26527544[accn])
- ZM1: [SRX26527541](https://www.ncbi.nlm.nih.gov/sra/SRX26527541[accn]) 
- ZM2: [SRX26527542](https://www.ncbi.nlm.nih.gov/sra/SRX26527542[accn])

### WGS BAM
CBM2 and CBM3 refer to two biological replicates. Each sequenced twice (1 and 2). 

- CBM2_1: [SRX25617961](https://www.ncbi.nlm.nih.gov/sra/SRX25617961[accn])
- CBM2_2: [SRX25617962](https://www.ncbi.nlm.nih.gov/sra/SRX25617962[accn])
- CBM3_1: [SRX25617963](https://www.ncbi.nlm.nih.gov/sra/SRX25617963[accn])
- CBM3_2: [SRX25617964](https://www.ncbi.nlm.nih.gov/sra/SRX25617964[accn])

### Nanopore hMeDIP-seq
Each numbered file represents a technical replicate from a single mouse (CBM1). 

- input: [SRX25617968](https://www.ncbi.nlm.nih.gov/sra/SRX25617968[accn])
- 1: [SRX25617965](https://www.ncbi.nlm.nih.gov/sra/SRX25617965[accn])
- 2: [SRX25617966](https://www.ncbi.nlm.nih.gov/sra/SRX25617966[accn])
- 3: [SRX25617967](https://www.ncbi.nlm.nih.gov/sra/SRX25617967[accn])

### If downloading bam format or aligning from basecalling
Pre-processing is required first by sorting according to chromosomal position (`samtools sort`) and indexing using `samtools index`.
 
### Extracting modified bases
Extracting modified bases was performed using `modkit pileup` and `modkit pileup-hemi` using the `--cpg` and `--mask` options, using the UCSC mm39 reference genome as a reference. 

Where read information is required `modkit extract calls --cpg --pass-only --reference mm39.fa` is used. This is specified within scripts. 

## Machine data
We recommend only downloading the processed files. These raw machine data are very large and take a long time to base-call and process. If needed however, they can be downloaded in fast5 format can be downloaded from here:

- CBM2_1: [SRX26304802](https://www.ncbi.nlm.nih.gov/sra/SRX26304802[accn])
- CBM2_2: [SRX26304803](https://www.ncbi.nlm.nih.gov/sra/SRX26304803[accn])
- CBM3_1: [SRX26304804](https://www.ncbi.nlm.nih.gov/sra/SRX26304804[accn])
- CBM3_2: [SRX26304805](https://www.ncbi.nlm.nih.gov/sra/SRX26304805[accn])

Machine data is available in pod5 format for the hMeDIP-seq replicates on Zenodo: https://doi.org/10.5281/zenodo.14514704. Please note that, due to individual file size limits, the archive for hMeDIP replicate 3 had to be broken into 5 smaller archives. Simply untar these into the same output directory and basecall using the `--recursive` option to include all data.   

__Raw machine data are not available for any of the Zymo DNA Methylation Standard samples.__

### If downloading raw fast5
Fast5 data requires conversion to pod5 format file using pod5 tools before base-calling and aligning (to the mm39/GRCm39 reference genome) with Dorado. Do so using the pod5-file-format package: https://github.com/nanoporetech/pod5-file-format. 