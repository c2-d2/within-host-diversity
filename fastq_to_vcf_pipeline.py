"""
Michael Martin mim3073@mail.harvard.edu

The purpose of this very simple script is to take input fastq files and do the following:


List of things:
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
-trimmomatic (.jar file needs to be in the same directory you are running the script from)
-bedtools
-bwa
-samtools
-gatk (.jar file needs to be in the same directory you are running the script from)
-picard (.jar file needs to be in the same directory you are running the script from)
-snpeff (.jar file and snpEff.config needs to be in teh same directory you are running the script from)

Furthermore, you need some files for the script to run as intended:
-.bed file of PE/PPE genes called peppe.bed (needs to be in the same directory you are running the script from)
-tab delimited file of drug resistance conferring mutations called dr_mutations.tab (gene \t mutation \t drug) (needs to be in the same directory you are running the script from)
-.fasta file to be used as a reference

When running the script, the first argument needs to be the reference .fasta, the second argument is the read frequency
cutoff for calling an hSNP (1 - threshold = cSNP cutoff). Remaining arguments are paired end.fastq reads.

While there is some flexibility in the reference used, some code elements are designed to specifically work with
H37RV. For example, when changing the chromosome name to match that in SnpEff, 'NC_000962.3' specifically refers to
the H37Rv reference. Also, the SnpEff database which is referenced is the H37Rv database. Both of these references should
be very easy to change for the reference of choice.

"""
#sysarg=script,ref,rare variant read frequency cutoff, samples
import sys
import os
import csv

# Stores the name of the reference fasta
reference=sys.argv[1]

# Stores the read frequency cutoff for consensus SNPs
snpcutoff=1-float(sys.argv[2])

# Stores the read frequency cutoff for heterogeneous SNPs
readcutoff=float(sys.argv[2])
print(sys.argv)
# Stores the names of all input samples
samples=[]
for i in xrange(3,len(sys.argv)):
	if (sys.argv[i].split('_')[0]) not in samples:
		samples.append(sys.argv[i].split('_')[0])

# 1. Prepares the reference genome to be used with BWA and GATK
# Creates BWA index
#os.system("bwa index "+reference)
# Generates fasta index file
#os.system('samtools faidx '+reference)

# Generates sequence dictionary
#os.system(
#	"java -jar picard.jar CreateSequenceDictionary R="+reference+" O="+reference.split('.fasta')[0]+".dict")
print(samples)
# Loops through the list of samples
trimmed_samples=[]

for sample in samples:
	print sample
	sample_trimmed=sample+'trimmed'
	trimmed_samples.append(sample_trimmed)
	'''
	# 2. Trims the raw fastq files
	os.system(
		'java -jar trimmomatic-0.36.jar PE -phred33 '+sample+'_1.fastq.gz '+sample+'_2.fastq.gz '+
		sample+'trimmed_1.fastq.gz '+sample+'trimmed_1_unpaired.fastq.gz '+sample+'trimmed_2.fastq.gz '+
		sample+'trimmed_2_unpaired.fastq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36')



	# 3. Aligns the reads to the reference genome and sorts the BAM
	os.system("bwa mem -M "+reference+" "+sample_trimmed+"_1.fastq.gz "+sample_trimmed+"_2.fastq.gz | samtools sort -O bam -o "+sample_trimmed+
			  "_bwa_sorted.bam")

	# 4. Add read groups to sorted BAM
	os.system("java -jar picard.jar AddOrReplaceReadGroups INPUT="+sample_trimmed+"_bwa_sorted.bam OUTPUT="+sample_trimmed+
			  "_bwa_sorted_RG.bam SORT_ORDER=coordinate RGID="+sample_trimmed+" RGLB=unknown RGPL=Illumina RGSM="+sample_trimmed+
			  " RGPU=project CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT")

	# 5. Identify target intervals for local realignment
	os.system("java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R "+reference+" -I "+sample_trimmed+
			  "_bwa_sorted_RG.bam -o "+sample_trimmed+"_bwa_sorted_RG_target.intervals")

	# 6 Conduct local realignment and filter out reads <100BP (soft and hard) on BAM and identify regions with no coverage
	os.system("java -jar GenomeAnalysisTK.jar -rf OverclippedRead --filter_is_too_short_value 100 --do_not_require_softclips_both_ends -rf ReadLength -minRead 100 -maxRead 500 -T IndelRealigner -R "+
			  reference+" -I "+sample_trimmed+"_bwa_sorted_RG.bam -targetIntervals "+sample_trimmed+"_bwa_sorted_RG_target.intervals -o "+
			  sample_trimmed+"_bwa_sorted_RG_RA_100.bam")
	os.system("bedtools genomecov -ibam "+sample_trimmed+"_bwa_sorted_RG_RA_100.bam  -bga | awk '$4==0' > "+
			  sample_trimmed+"_cov.bed")

	# Identifies regions with no coverage
	os.system("bedtools genomecov -ibam "+sample_trimmed+"_bwa_sorted_RG_RA_100.bam  -bga | awk '$4==0' > "+
			  sample_trimmed+"_cov_100.bed")

	# 7. Create mpileup from BAM
	os.system("samtools mpileup -d 1000 -q 30 -t DP,AD,ADF,ADR,SP -u -g -f "+reference+" "+sample_trimmed+
			  "_bwa_sorted_RG_RA_100.bam -o "+sample_trimmed+"_bwa_sorted_RG_RA_100_mq30_baq.mpileup")

	# 8. Call and filter variants from mpileup
	os.system("bcftools call -O v -m -v "+sample_trimmed+
			  "_bwa_sorted_RG_RA_100_mq30_baq.mpileup | bcftools filter -i 'SP<60 & ADF[1]>1 & ADR[1]>1 & MQ>30 & QUAL>50 & FORMAT/DP > 20 & SP<60 & ADF[1]>1 & ADR[1] >1' -o "+
			  sample_trimmed+"_bwa_sorted_RG_RA_100_mq30_baq_bq50_dp20_sp60_ad1.vcf")

	# Filters any regions w/in 12 bp of two cSNPs
	# cSNP VCF
	os.system("bcftools call -O v -m -v " + sample_trimmed +
			  "_bwa_sorted_RG_RA_100_mq30_baq.mpileup | bcftools filter -i '(AD[1]/(AD[0]+AD[1]) >= "+str(snpcutoff)+") & SP<60 & ADF[1]>1 & ADR[1]>1 & MQ>30 & QUAL>50 & FORMAT/DP > 20' -o " +
			  sample_trimmed + "_bwa_sorted_RG_RA_100_mq30_baq_bq50_dp20_sp60_ad1_cons.vcf")
	csnps = list(csv.reader(open(sample_trimmed + '_bwa_sorted_RG_RA_100_mq30_baq_bq50_dp20_sp60_ad1_cons.vcf', 'rb'),
							delimiter='\t'))
	# Finds where the header ends
	for k in xrange(0, 100):
		if csnps[k][0] == '#CHROM':
			startindex = k
			break

	csnp_pos = []
	# For each row after the header
	for j in xrange(startindex + 1, len(csnps)):
		# if this position is a SNP (one reference base and one variant base)
		if len(csnps[j][3]) == 1 and len(csnps[j][4]) == 1:
			csnp_pos.append(int(csnps[j][1]))
	excluded_regions = []
	for j in xrange(0, len(csnp_pos) - 1):
		if csnp_pos[j + 1] - csnp_pos[j] <= 12:

			# Region between the two variants
			# Ensures size of the region is positive
			if csnp_pos[j + 1] - csnp_pos[j] > 1:
				excluded_regions.append(['NC_000962.3'])
				excluded_regions[-1].append(csnp_pos[j] + 1)
				excluded_regions[-1].append(csnp_pos[j + 1] - 1)
				excluded_regions[-1].append('0')

			# Region before the first variant which is still w/in 12BP of second variant
			# Ensures size of the region is positive
			if (csnp_pos[j + 1] - 12) < (csnp_pos[j] - 1):
				excluded_regions.append(['NC_000962.3'])
				excluded_regions[-1].append(csnp_pos[j + 1] - 12)
				excluded_regions[-1].append(csnp_pos[j] - 1)
				excluded_regions[-1].append('0')

			# Region after the second variant which is still w/in 12BP of first variant
			# Ensures size of the region is positive
			if (csnp_pos[j + 1] + 1) < (csnp_pos[j] + 12):
				excluded_regions.append(['NC_000962.3'])
				excluded_regions[-1].append(csnp_pos[j + 1] + 1)
				excluded_regions[-1].append(csnp_pos[j] + 12)
				excluded_regions[-1].append('0')

	with open(sample_trimmed + "_density_filtered.bed", 'w') as file:
		writer = csv.writer(file, delimiter='\t')
		writer.writerows(excluded_regions)

	# Removes any SNPs w/in 12 bp of two cSNPs
	os.system("bedtools intersect -v -header -a " +
			  sample_trimmed + "_bwa_sorted_RG_RA_100_mq30_baq_bq50_dp20_sp60_ad1.vcf -b " +
			  sample_trimmed + "_density_filtered.bed -wa > " + sample_trimmed + "_bwa_sorted_RG_RA_100_mq30_baq_bq50_dp20_sp60_ad1_dens.vcf")

	# Finds positions which fail the SP, MQ, DP, and Qual thresholds
	os.system("bcftools call -O v -m -v " + sample_trimmed +
			  "_bwa_sorted_RG_RA_100_mq30_baq.mpileup | bcftools filter -i 'SP>=60 || MQ<=30 || FORMAT/DP<=20 || QUAL<=50' -o " +
			  sample_trimmed + "_bwa_sorted_RG_RA_100_mq30_baq_bq50_dp20_sp60_ad1_failed_filter.vcf")
	# Finds positions which fail the AD threshold
	os.system("bcftools call -O v -m -v " + sample_trimmed +
			  "_bwa_sorted_RG_RA_100_mq30_baq.mpileup | bcftools filter -i '(AD[1]/(AD[0]+AD[1]) > "+str(readcutoff)+") & (ADF[1]<=1 || ADR[1]<=1)' -o " +
			  sample_trimmed + "_bwa_sorted_RG_RA_100_mq30_baq_bq50_dp20_sp60_ad1_failed_ad_filter.vcf")
'''
# 9. Removes variants which fail quality filters in any samples, have coverage of 0 in any samples, and which are in PE/PPE regions
# Makes a .bed file called 'filtered_sites.bed' which contains the positions of variants which failed quality filters in
# any of the samples.
nocov=[]
sample_trimmedbeds=''
excl_poss=[]
excl_pos=[]
poss=[]
for sample in samples:
	sample_trimmedbeds+=(sample_trimmed+"_cov_100.bed ")

	list1 = list(csv.reader(open(sample_trimmed+'_bwa_sorted_RG_RA_100_mq30_baq_bq50_dp20_sp60_ad1_failed_filter.vcf', 'rb'), delimiter='\t'))
	# Finds where the header ends
	for k in xrange(0, 100):
		if list1[k][0]=='#CHROM':
			startindex = k
			break
	# For each row after the header
	for j in xrange(startindex + 1, len(list1)):
		# If this position is a SNP (one reference base and one variant base)
		if len(list1[j][3]) == 1 and len(list1[j][4]) == 1:
			if list1[j][1] not in poss:
				poss.append(list1[j][1])
				excl_pos.append('NC_000962.3')
				excl_pos.append(int(list1[j][1]))
				excl_pos.append(int(list1[j][1]))
				excl_pos.append('0')
				excl_poss.append(excl_pos)
				excl_pos = []
	list1 = list(csv.reader(open(sample_trimmed+'_bwa_sorted_RG_RA_100_mq30_baq_bq50_dp20_sp60_ad1_failed_ad_filter.vcf', 'rb'), delimiter='\t'))
	# Finds where the header ends
	for k in xrange(0, 100):
		if '#CHROM' in list1[k]:
			startindex = k
			break
	# For each row after the header
	for j in xrange(startindex + 1, len(list1)):
		# If this position is a SNP (one reference base and one variant base)
		if len(list1[j][3]) == 1 and len(list1[j][4]) == 1:
			if list1[j][1] not in poss:
				poss.append(list1[j][1])
				excl_pos.append('NC_000962.3')
				excl_pos.append(int(list1[j][1]))
				excl_pos.append(int(list1[j][1]))
				excl_pos.append('0')
				excl_poss.append(excl_pos)
				excl_pos = []

with open("filtered_sites.bed", 'w') as file:
	writer = csv.writer(file, delimiter='\t')
	writer.writerows(excl_poss)

# Combines the sites which failed filters into the same .bed file as the PE/PPE regions and regions with no coverage
os.system("cat filtered_sites.bed peppe.bed "+sample_trimmedbeds+" > remove_sites.bed")

for sample_trimmed in trimmed_samples:
	# Removes positions which failed filters
	os.system("bedtools intersect -v -header -a "+sample_trimmed+"_bwa_sorted_RG_RA_100_mq30_baq_bq50_dp20_sp60_ad1_dens.vcf -b remove_sites.bed -wa > "+sample_trimmed+"_bwa_sorted_RG_RA_100_mq30_baq_bq50_dp20_sp60_ad1_cov_pe_dens_passfilter.vcf")

	# 10. Annotates the final VCF
	# Changes chromosome name in VCF to match that in SnpEff database
	os.system("sed 's/NC_000962.3/Chromosome/' "+sample_trimmed+"_bwa_sorted_RG_RA_100_mq30_baq_bq50_dp20_sp60_ad1_cov_pe_dens_passfilter.vcf >"+sample_trimmed+"_bwa_sorted_RG_RA_100_mq30_baq_bq50_dp20_sp60_ad1_cov_pe_dens_passfilter_chr.vcf")
	os.system('java -jar SnpEff.jar Mycobacterium_tuberculosis_h37rv '+sample_trimmed+'_bwa_sorted_RG_RA_100_mq30_baq_bq50_dp20_sp60_ad1_cov_pe_dens_passfilter_chr.vcf >'+sample_trimmed+'_bwa_sorted_RG_RA_100_mq30_baq_bq50_dp20_sp60_ad1_cov_pe_dens_passfilter_ann.vcf')
	# Removes the VCF with the chromosome changed
	os.system('rm '+sample_trimmed+'_bwa_sorted_RG_RA_100_mq30_baq_bq50_dp20_sp60_ad1_cov_pe_dens_passfilter_chr.vcf')

# Opens the tab delimited list of mutations
mutations=list(csv.reader(open('dr_mutations.tab','rb'),delimiter='\t'))
this_gene=[]
dr_mutations=[]
this_mutation=[]

# Identifies mutations in the VCFs which are associated with resistance
for sample_trimmed in trimmed_samples:
    list1=list(csv.reader(open(sample_trimmed+'_bwa_sorted_RG_RA_100_mq30_baq_bq50_dp20_sp60_ad1_cov_pe_dens_passfilter_ann.vcf','rb'),delimiter='\t'))
    if sample_trimmed not in dr_mutations:
        dr_mutations.append(sample_trimmed)
    for k in xrange(0,len(list1)):
        for j in xrange(0,len(mutations)):
            if any(mutations[j][0] in s for s in list1[k]):
                this_gene.append(list1[k])
                for l in xrange(0,len(this_gene)):
                    if any(mutations[j][1] in s for s in this_gene[l]):
                        this_mutation.append(mutations[j][0])
                        this_mutation.append(list1[k][3]+list1[k][1]+list1[k][4])
                        this_mutation.append(mutations[j][1])
                        this_mutation.append(mutations[j][2])
                        splitted=list1[k][-1].split(':')
            			# Splits the second to last item in the list (FORMAT), this is where
            			# the indices for the data comes from
            			# GT:PL:DP:SP:ADF:ADR:AD
                        splitted2=list1[k][-2].split(':')
            			# Supporting variant read frequency
                        this_mutation.append(str(float(splitted[splitted2.index('AD')].split(',')[1])/float(splitted[splitted2.index('DP')])))
                        dr_mutations.append(this_mutation)
                        this_mutation=[]
                this_gene=[]
    dr_mutations.append('\n')

with open('detected_dr_mutations.tab', 'w') as file:
    file.writelines('\t'.join(i) + '\n' for i in dr_mutations)
