# Within Host Diversity

These are a set of simple scripts designed to call variants from short read whole genome sequencing M. tuberuclosis specimens. The emphasis is on the identification of loci with heterogeneous base calls with the goal of identifying within-host diversity to be used in transmision estimates.

There are two primary scripts, both written in Python v2.7: 
1. fastq_to_vcf_pipeline.py
	This script is designed to take as an input raw paired end illumina reads and output filtered VCF files for variant positions as well as a tab delimited file with the identified drug resistant mutations in each sample. 

	The script performs the following steps:
	1. Prepares the reference for use with BWA and GATK
	2. Trims the reads with trimmomatic
	3. Aligns the reads to the reference genome and sorts the BAM
	4. Add read groups to the sorted BAM
	5. Identify target intervals for local realignment
	6. Conduct local realignment and filter out reads <100 bases and identify regions with no coverage
	7. Create mpileup from BAM
	8. Call and filter variants from mpileup
	9. Remove variants which fail quality thresholds in any samples and which occur in PE/PPE regions
	10. Identify drug resistance conferring mutations

	In order to run properly, this script requires a number of bioinformatics tools:
	* Trimmomatic (.jar file needs to be in the same directory you are running the script from) (http://www.usadellab.org/cms/?page=trimmomatic)
	* BEDtools (http://bedtools.readthedocs.io/en/latest/)
	* BWA (http://bio-bwa.sourceforge.net)
	* SAMtools/BCFtools (http://www.htslib.org)
	* GATK3 (.jar file needs to be in the same directory you are running the script from) (https://gatkforums.broadinstitute.org/gatk)
	* Picard (.jar file needs to be in the same directory you are running the script from) (https://broadinstitute.github.io/picard/index.html)
	* SnpEff (http://snpeff.sourceforge.net)

	Furthermore, you need some files for the script to run as intended:
	* .bed file of PE/PPE genes called peppe.bed (needs to be in the same directory you are running the script from) (Camus et al., 2002)
	* Tab delimited file of drug resistance conferring mutations called dr_mutations.tab (gene \t mutation \t drug) (needs to be in the same directory you are running the script from) (Table S8.1, Walker et al. 2015)
	* Fasta file to be used as a reference (GenBank accession number: NC_000962.3)


2. within_host_heterogeneity.py
	This script is designed to take as input variant position VCF files and use them to pull out variant sites with a focus on finding sites with heterogeneous variant calls. This script only takes SNPs into account and will remove all indels. 

	The script will produce a number of outputs: 

    * cSNP and hSNP alignment
    * Pairwise cSNP distances
    * Pairwise cSNP loci
    * List of SNPs, hSNPs, cSNPs
    * SNP distance from the reference
    * Informative hSNP loci with read frequency support
    * Plot of pairwise cSNP distances
    * Histogram of stand bias scores
    * Histogram of hSNP shared between plot
    * Plot with DP, BQ, and MQ
    * Will also output information to std out:
        -Correlation between cSNPs from reference and number of hSNPs
        -Minimum, maximum, mean, and std values comparing cSNPs to hSNPs for DP, QUAL, MQ, SP, BQB, MQB, RPB
        -Maximum distance from reference genome
        -Minimum, maximum, mean, std pairwise SNP distances

cSNP = consensus SNP above snpcutoff read read_frequency
hSNP = heterogeneous SNP where the variant read frequency is between the
use defined threshold and the snpcutoff

When running the script, the first argument should be the threshold for calling hSNPs (in other words, 0.10 would mean
that any variant w/ <=10% supporting reads is called as reference and any variant with >10% and >=90% supporting reads
is called as a consensus SNP.

## Reference
>*In preperation*

## Contact
>mim3073@mail.harvard.edu