"""
Michael Martin mim3073@mail.harvard.edu

The purpose of this script is to take input VCF files and use them to pull out variant sites with a focus on finding
sites with heterogeneous variant calls. This script only takes SNPs into account and will remove all indels.

The script will calculate and output a number of files:

-cSNP and hSNP alignment
-Pairwise cSNP distances
-Pairwise cSNP loci
-List of SNPs, hSNPs, cSNPs
-SNP distance from the reference
-Informative hSNP loci with read frequency support
-Plot of pairwise cSNP distances
-Histogram of stand bias scores
-Histogram of hSNP shared between plot
-Plot with DP, BQ, and MQ
-Will also output information to std out:
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

"""


# import required modules
import csv
import sys
from sys import argv
import matplotlib.pyplot as plt
from datetime import datetime
import seaborn as sns
from scipy import stats
import numpy as np
import pandas as pd

# Times how long the script takes to run
startTime = datetime.now()

# Stores the cutoff for calling variants as reference and the cutoff for calling a cSNP
readcutoff=float(sys.argv[1])
snpcutoff=(1-float(sys.argv[1]))


# define required lists
res=[]
snps=[]
matrixed=[]
splitted=[]
combo=[]
samples=[]
seqs=[]
missingseq=[]
distoutputmatrix=[[] for x in xrange(3,len(sys.argv))]
idoutputmatrix=[[] for x in xrange(3,len(sys.argv))]
csnp_dist=[]
csnp_dp=[]
hsnp_dp=[]
csnp_sp=[]
hsnp_sp=[]
csnp_qual=[]
hsnp_qual=[]
csnp_rpb=[]
hsnp_rpb=[]
csnp_mq=[]
hsnp_mq=[]
csnp_bqb=[]
hsnp_bqb=[]
csnp_mqb=[]
hsnp_mqb=[]
csnps=[]
hsnps=[]
snpindex=[]
this_hsnp=[]
this_sample_this_hsnp=[]
this_sample_this_hsnp_cont=[]
temp_hsnps=[]
hsnp_cont=[]


# Read in the VCF input, excludes any variants at a read frequency <threshold (set by arg 1)
# Data for each variant is stored in for order: read frequency, DP, QUAL, MQ, SP, BQB, MQB, RPB
for i in xrange(3,len(sys.argv)):
	# Stores list of samples
	samples.append(sys.argv[i][0:6])
	# Opens VCF
	list1=list(csv.reader(open(sys.argv[i],'rb'),delimiter='\t'))
	# Finds where the header ends
	for k in xrange(0,100):
		if '#CHROM' in list1[k]:
			startindex=k
	# for each row after the header
	for j in xrange(startindex+1,len(list1)):
		# if this position is a SNP (one reference base and one variant base) as opposed ot an indel
		if len(list1[j][3])==1 and len(list1[j][4])==1:
			# Defines needed lists
			data=[]
			reads=[]
			splitted=[]
			splitted2=[]
			# splits the last item in the list (FORMAT), this is where the data is:
			splitted=list1[j][len(list1[j])-1].split(':')
			# Splits the second to last item in the list (FORMAT), this is where
			# the indices for the data comes from
			# GT:PL:DP:SP:ADF:ADR:AD
			splitted2=list1[j][len(list1[j])-2].split(':')
			# Splits the INFO list
			splitted3=list1[j][7].split(';')
			splitted3_label=list(item.split('=')[0] for item in splitted3)
			splitted3_data=list(item.split('=')[1] for item in splitted3)
			# Position of the variant
			reads.append(int(list1[j][1]))
			# Variant in form ref+position+alt
			reads.append(list1[j][3]+list1[j][1]+list1[j][4])
			# ample name
			reads.insert(2,sys.argv[i][0:6])
			# Data in the format: RF, DP, QUAL, MQ, SP, BQB, MQB, RPB
			# Read Frequency
			data.append(float(splitted[splitted2.index('AD')].split(',')[1])/\
			float(splitted[splitted2.index('DP')]))
			# DP
			data.append(float(splitted[splitted2.index('DP')]))
			# QUAL
			data.append(float(list1[j][5]))
			# MQ
			data.append(float(splitted3_data[splitted3_label.index('MQ')]))
			# SP
			data.append(float(splitted[splitted2.index('SP')]))

			# BQB
			if 'BQB' in splitted3_label:
				data.append(float(splitted3_data[splitted3_label.index('BQB')]))
			# If 100% variant reads, assume BQB=1
			else:
				data.append(1.0)

			# MQB
			if 'MQB' in splitted3_label:
				data.append(float(splitted3_data[splitted3_label.index('MQB')]))
			# If 100% variant reads, assume MQB=1
			else:
				data.append(1.0)

			# RPB
			if 'RPB' in splitted3_label:
				data.append(float(splitted3_data[splitted3_label.index('RPB')]))
			# If 100% variant reads, assume RPB=1
			else:
				data.append(1.0)
			reads.append(data)
			# If this position has the necessary supporting variant read frequency to be called an hSNP/cSNP
			if reads[-1][0]>readcutoff:
				combo.append(reads)



# Combine items in list that are the same variant
for row in combo:
	for i in xrange(0,len(res)):
		if row[1]==res[i][1]:
			res[i] += row[2:]				
			break
	else:
		res.append(row)

# Adds spacers for column alignment in output csv
for i in xrange(0,len(res)):
	for j in xrange(2,len(res[i]),2):
		seqs.append(res[i][j])
	for l in xrange(0,len(samples)):
		if samples[l] not in seqs:
			missingseq.append(samples[l])
	for m in xrange(0,len(missingseq)):
		for k in xrange(0,2):
			res[i].insert(2+(samples.index(missingseq[m]))*2,"")
	snps.append(res[i])
	missingseq=[]
	seqs=[]

# Sorts the list by the position of the variant
snps.sort(key=lambda x: int(x[0]))

# Iterates through the list of SNPs to output a file with the informative hSNPs and their supporting read frequency by
# sample and to output list of cSNPs and hSNPs. Also tabulates cSNP and hSNP statistics.
this_consensus_snp=[]
hsnp_shared_between=[]
for i in xrange(0,len(snps)):
	heterogeneous_site_counter=0
	consensus_site_counter=0
	this_hsnp_pos=[]
	temp_hsnps.append(snps[i][0])
	temp_hsnps.append(snps[i][1])
	this_consensus_snp.append(snps[i][0])
	this_consensus_snp.append(snps[i][1])

	# Iterates through each individual snp list and appends the read frequency
	# from each sample at this position to read_frequencies

	for k in xrange(3,len(snps[i]),2):
		# If there is a variant for this sample at this position
		if len(snps[i][k])>0:
			# If the variant is less than the threshold for calling a cSNP
			if snps[i][k][0]<snpcutoff:

				# This is for the hsnp_bases_cont output file
				# Adds the supporting variant read frequency for this sample
				this_sample_this_hsnp_cont.append(float(snps[i][k][0])*100)

				# These are for the hsnp_list output file and for creating the hSNP alignments
				# Appending the sample name as well as the list of data
				temp_hsnps.append(snps[i][k-1])
				temp_hsnps.append(snps[i][k])

				# These are for the csnp list output file and for creating the cSNP alignments
				# Blank spacers because no cSNP here
				this_consensus_snp.append('')
				this_consensus_snp.append('')

				# hSNP statistics
				hsnp_dp.append(snps[i][k][1])
				hsnp_qual.append(snps[i][k][2])
				hsnp_mq.append(snps[i][k][3])
				hsnp_sp.append(snps[i][k][4])
				hsnp_bqb.append(snps[i][k][5])
				hsnp_mqb.append(snps[i][k][6])
				hsnp_rpb.append(snps[i][k][7])

				#Tabulates the number of hSNPs at this loci
				heterogeneous_site_counter+=1

			# If the variant is more than the threshold for calling a cSNP
			elif snps[i][k][0]>=snpcutoff:


				# This is for the hSNP_bases_cont output file
				# Adds the supporting variant read frequency
				this_sample_this_hsnp_cont.append(float(snps[i][k][0])*100)

				# These are for the hSNP_list output file and for creating the hSNP alignments
				# Adds the sample name as well as the list of data
				temp_hsnps.append(snps[i][k-1])
				temp_hsnps.append(snps[i][k])

				# These are for the cSNP list output file and for creating the cSNP alignments
				# Adds the sample name as well as th elist of data
				this_consensus_snp.append(snps[i][k-1])
				this_consensus_snp.append(snps[i][k])

				#cSNP statistics
				csnp_dp.append(snps[i][k][1])
				csnp_qual.append(snps[i][k][2])
				csnp_mq.append(snps[i][k][3])
				csnp_sp.append(snps[i][k][4])
				csnp_bqb.append(snps[i][k][5])
				csnp_mqb.append(snps[i][k][6])
				csnp_rpb.append(snps[i][k][7])

				# Tabulates the number of cSNPs at this loci
				consensus_site_counter+=1

		# If there is not a variant for this sample at this position
		else:

			# These are for the hSNP_list output file and for creating the hSNP alignments
			# Blank spacers because no hSNP here
			temp_hsnps.append('')
			temp_hsnps.append('')

			# These are for the cSNP list output file and for creating the cSNP alignments
			# Blank spacers because no cSNP here
			this_consensus_snp.append('')
			this_consensus_snp.append('')

			#This is for the hSNP_bases_cont output file
			this_sample_this_hsnp_cont.append(0.0)

		# Adds the list for this sample to the list for this variant


	# Makes sure that there is at least one heterogeneous variant at this loci
	if 0<heterogeneous_site_counter:
		# These are for the hSNP_list output file and for creating the hSNP alignments
		hsnps.append(temp_hsnps)

		# Calculates how many samples this loci is an hSNP or a cSNP in
		hsnp_shared_between.append(len(samples)-this_sample_this_hsnp_cont.count(0.0))

		# Makes sure that this loci isn't variant in every sample
		if heterogeneous_site_counter<len(samples):

			# This is for the hSNP_bases_cont output file
			# Add the position of the variant
			this_sample_this_hsnp_cont.insert(0, snps[i][1])
			# Adds this position to the master lists
			hsnp_cont.append(this_sample_this_hsnp_cont)

	# This is for the cSNP_list output file
	# Makes sure there is at least one consensus variant at this site
	if 0<consensus_site_counter:
		csnps.append(this_consensus_snp)

	this_hsnp.append(this_sample_this_hsnp)
	this_consensus_snp=[]
	temp_hsnps=[]
	this_sample_this_hsnp_cont=[]
	this_hsnp=[]

#Adds the list of samples to the hsnp_bases_cont output file
hsnp_cont.insert(0,[''])
hsnp_cont[0].extend(samples)

# IUPAC bases
R=['A','G']
Y=['C','T']
W=['A','T']
S=['G','C']
M=['A','C']
K=['G','T']

hsnp_counter=[]
alignment=[]
thisalignment=[]
this_hsnp_alignment=[]
hsnp_alignment=[]
num_snps=[]
num_csnps=[]
num_hsnps=[]
snp_dist=[]

# Calculates pairwise cSNP distances and makes cSNP and hSNP alignments
# Iterates through the list of samples
for l in xrange(0,len(samples)):
	hsnp_counter = 0
	csnp_counter = 0
	# Iterates through the list of samples
	for j in xrange(0,len(samples)):
		# Iterates through the list of SNPs
		for i in xrange(0,len(csnps)):
			# If only one of the samples has a variant, ensure it is >=snpcutoff read frequency
			if (samples[l] in csnps[i]) and (samples[j] not in csnps[i]) and \
					(float(csnps[i][csnps[i].index(samples[l])+1][0])>=snpcutoff):
				matrixed.append(str(csnps[i][1])+samples[l]+samples[j])
				if snps[i][0] not in snpindex:
					snpindex.append(csnps[i][0])
			if (samples[j] in csnps[i]) and (samples[l] not in csnps[i]) and \
					(float(csnps[i][csnps[i].index(samples[j])+1][0])>=snpcutoff):
				matrixed.append(str(csnps[i][1])+samples[l]+samples[j])
				if csnps[i][0] not in snpindex:
					snpindex.append(snps[i][0])

			# If both samples have a snp, but only one is >=snpcutoff read frequency
			if (samples[l] in csnps[i]) and (samples[j] in csnps[i]) and \
					(float(csnps[i][csnps[i].index(samples[l])+1][0])>=snpcutoff) and \
					(float(csnps[i][csnps[i].index(samples[j])+1][0])<snpcutoff):
				matrixed.append(str(csnps[i][1])+samples[l]+samples[j])
				if csnps[i][0] not in snpindex:
					snpindex.append(snps[i][0])
			if (samples[l] in csnps[i]) and (samples[j] in csnps[i]) and \
					(float(csnps[i][csnps[i].index(samples[l])+1][0])<snpcutoff) and \
					(float(csnps[i][csnps[i].index(samples[j])+1][0])>=snpcutoff):
				matrixed.append(str(snps[i][1])+samples[l]+samples[j])
				if csnps[i][0] not in snpindex:
					snpindex.append(csnps[i][0])

		distoutputmatrix[l].append(len(matrixed))

		if l>j:
			csnp_dist.append(len(matrixed))
		idoutputmatrix[l].append((snpindex))

		snpindex=[]
		matrixed=[]


	# cSNP and hSNP alignments
	# For each sample add a fasta header
	alignment.append(">" + samples[l][0:6])
	hsnp_alignment.append(">" + samples[l][0:6])
	# Iterates through list of cSNPs
	for j in xrange(0, len(csnps)):
		# If this snp is present in this sample
		if samples[l] in csnps[j]:
			# If the read frequency is >=snpcutoff add the variant base to the alignment
			thisalignment.extend(snps[j][1][-1])
			csnp_counter+=1
		# Else add the reference base to the alignment
		else:
			thisalignment.extend(csnps[j][1][0])

	# Combine the alignment to the master list
	alignment.append(''.join(map(str, thisalignment)))
	thisalignment = []

	# Iterates through list of hSNPs
	for j in xrange(0, len(hsnps)):

		# If this SNP is present in this sample
		if samples[l] in hsnps[j]:
			hsnp_counter+=1

			these_bases = []
			these_bases.append(hsnps[j][1][0])
			these_bases.append(hsnps[j][1][-1])
			# If this SNP is an hSNP
			if hsnps[j][hsnps[j].index(samples[l]) + 1][0] < snpcutoff:

				# Figure out which ambiguous base should be used
				if len(set(these_bases).intersection(R)) == 2:
					this_hsnp_alignment += 'R'
				if len(set(these_bases).intersection(Y)) == 2:
					this_hsnp_alignment += 'Y'
				if len(set(these_bases).intersection(W)) == 2:
					this_hsnp_alignment += 'W'
				if len(set(these_bases).intersection(S)) == 2:
					this_hsnp_alignment += 'S'
				if len(set(these_bases).intersection(M)) == 2:
					this_hsnp_alignment += 'M'
				if len(set(these_bases).intersection(K)) == 2:
					this_hsnp_alignment += 'K'

			# If this SNP is a cSNP add just the variant base
			elif hsnps[j][hsnps[j].index(samples[l]) + 1][0] >= snpcutoff:
				this_hsnp_alignment += hsnps[j][1][-1]
				csnp_counter += 1
		# Else add the reference
		else:
			this_hsnp_alignment += hsnps[j][1][0]

	# Combine the alignment to the master list
	hsnp_alignment.append(''.join(map(str, this_hsnp_alignment)))
	this_hsnp_alignment = []

	# Adds the number of cSNPs and hSNPs to a list
	sample_num_snps=[]
	num_hsnps.append(hsnp_counter)
	num_csnps.append(csnp_counter)
	sample_num_snps.append(samples[l])
	sample_num_snps.append(csnp_counter)
	sample_num_snps.append(hsnp_counter)
	num_snps.append(sample_num_snps)


# Adds row and column headers to the distance matrix
for i in xrange(0,len(distoutputmatrix)):
	distoutputmatrix[i].insert(0,samples[i])
distoutputmatrix.insert(0,[''])
distoutputmatrix[0].extend(samples)

# Adds referenceto alignment as outgroup
alignment.append(">REF")
for j in xrange(0,len(csnps)):
	thisalignment.append(csnps[j][1][0])
alignment.append(''.join(map(str,thisalignment)))

hsnp_alignment.append(">REF")
for j in xrange(0,len(hsnps)):
	this_hsnp_alignment.append(snps[j][1][0])
hsnp_alignment.append(''.join(map(str,this_hsnp_alignment)))

# Save cSNP and hSNP alignment
with open (argv[2]+"_cSNP_Alignment.aln",'w') as file:
	for i in xrange(0,len(alignment)-1):
		file.write(str(alignment[i])+'\n')
	file.write(str(alignment[len(alignment)-1]))

with open (argv[2]+"_hSNP_Alignment.aln",'w') as file:
	for i in xrange(0,len(hsnp_alignment)-1):
		file.write(str(hsnp_alignment[i])+'\n')
	file.write(str(hsnp_alignment[len(hsnp_alignment)-1]))

# Sorts the number of SNPs list by the cSNP distance to the references
num_snps.sort(key=lambda x: int(x[1]))

# Adds column headers
num_snps.insert(0,['Sample','cSNPs','hSNPs'])

for i in xrange(0,len(idoutputmatrix)):
	idoutputmatrix[i].insert(0,samples[i])
idoutputmatrix.insert(0,[''])
idoutputmatrix[0].extend(samples)

# Pairwise cSNP distances
with open (argv[2][0:len(argv[2])]+"_cSNP_Distance.csv",'w') as file:
    writer=csv.writer(file)
    writer.writerows(distoutputmatrix)

# cSNP ID
with open (argv[2][0:len(argv[2])]+"_cSNP_ID.csv",'w') as file:
    writer=csv.writer(file)
    writer.writerows(idoutputmatrix)

# SNP List
with open (argv[2]+"_SNP_List.csv",'w') as file:
	writer=csv.writer(file)
	writer.writerows(snps)

# cSNP List
with open (argv[2]+"_cSNP_List.csv",'w') as file:
	writer=csv.writer(file)
	writer.writerows(csnps)

# hSNP List
with open (argv[2]+"_hSNP_List.csv",'w') as file:
	writer=csv.writer(file)
	writer.writerows(hsnps)

# SNP distance from reference
with open (argv[2]+"_SNP_num.csv",'w') as file:
	writer=csv.writer(file)
	writer.writerows(num_snps)

# List of informative hSNP loci with read frequency support
with open (argv[2]+"_Inf_hSNP_loci_cont.csv",'w') as file:
	writer=csv.writer(file)
	writer.writerows([list(i) for i in zip(*hsnp_cont)])


#Makes Figures
#pretty colors
csnp_col = '#703D6F'
hsnp_col = '#6F8B94'
ref_hsnp_col='#706482'
my_grey = '#414451'
blue='#3E74CA'

# cSNP distance
fig=plt.figure(figsize=(12,9))
ax=fig.add_subplot(111)
ax.hist(csnp_dist,bins=100,normed=True,alpha=0.75,facecolor=csnp_col,label="cSNP")
ax.set_xlabel("cSNP Distance",fontsize=22,color=my_grey)
ax.set_ylabel("Frequency",fontsize=22,color=my_grey)
ax.tick_params(colors=my_grey,labelsize=18)
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color(my_grey)
ax.spines['left'].set_color(my_grey)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
plt.savefig(sys.argv[2]+"_cSNP_Distance",dpi=300)
plt.close()

# Histogram of strand bias scores
fig=plt.figure(figsize=(12,9))
ax=fig.add_subplot(111)
ax.hist(csnp_sp,normed=True,alpha=0.75,facecolor=csnp_col,label="cSNP")
ax.hist(hsnp_sp,normed=True,alpha=0.75,facecolor=hsnp_col,label="hSNP")
ax.set_xlabel("SP",fontsize=18,color=my_grey)
ax.set_ylabel("Frequency",fontsize=18,color=my_grey)
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color(my_grey)
ax.spines['left'].set_color(my_grey)
leg = ax.legend(loc="upper right", fontsize=16)
for text in leg.get_texts():
	plt.setp(text, color=my_grey)
plt.savefig(sys.argv[2]+"_SP",dpi=300)
plt.close()

# hSNP shared between histogram
fig=plt.figure(figsize=(12,9))
ax=fig.add_subplot(111)
ax.hist(hsnp_shared_between,normed=False,alpha=0.75,facecolor=my_grey,bins=25)
ax.set_xlabel("Number of Specimens",fontsize=18,color=my_grey)
ax.set_ylabel("Number of Shared hSNPs",fontsize=18,color=my_grey)
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.spines['bottom'].set_color(my_grey)
ax.spines['left'].set_color(my_grey)
plt.savefig(sys.argv[2]+"_hsnp_shared_between",dpi=300)
plt.close()

# DP, BQ, MQ Plot
bq_type=[]
bq=[]
for i in xrange(0,len(csnp_qual)):
	bq_type.append('cSNPs')
	bq.append(float(csnp_qual[i]))
for i in xrange(0,len(hsnp_qual)):
	bq_type.append('hSNPs')
	bq.append(float(hsnp_qual[i]))

bq_df = pd.DataFrame(
    {'Region': bq_type,
     'Base Quality': bq
    })


fig=plt.figure(figsize=(15,9))
ax1=fig.add_subplot(131)
ax1.hist(csnp_dp,normed=True,alpha=0.75,facecolor=csnp_col,label="cSNP",bins=25)
ax1.hist(hsnp_dp,normed=True,alpha=0.75,facecolor=hsnp_col,label="hSNP",bins=25)
ax1.set_xlabel("High Quality Read Depth",fontsize=18,color=my_grey)
ax1.set_ylabel("Frequency",fontsize=18,color=my_grey)
ax1.spines['top'].set_color('none')
ax1.spines['right'].set_color('none')
ax1.spines['bottom'].set_color(my_grey)
ax1.spines['left'].set_color(my_grey)
ax1.tick_params(labelcolor=my_grey,color=my_grey)
ax2=fig.add_subplot(132)
sns.stripplot(x="Region", y="Base Quality", data=bq_df,jitter=1,palette=[hsnp_col,csnp_col],order=['cSNPs','hSNPs'])
ax2.set_xlabel("Variant",fontsize=18,color=my_grey)
ax2.set_ylabel("Base Quality",fontsize=18,color=my_grey)
ax2.tick_params(labelcolor=my_grey)
ax2.set_ylim(0,250)
ax2.yaxis.set_ticks(np.arange(0, 251, 50))
ax2.spines['top'].set_color('none')
ax2.spines['right'].set_color('none')
ax2.spines['bottom'].set_color(my_grey)
ax2.spines['left'].set_color(my_grey)
ax3=fig.add_subplot(133)
ax3.hist(csnp_mq,normed=True,alpha=0.75,facecolor=csnp_col,label="cSNP",bins=25)
ax3.hist(hsnp_mq,normed=True,alpha=0.75,facecolor=hsnp_col,label="hSNP",bins=25)
ax3.set_xlabel("Mapping Quality",fontsize=18,color=my_grey)
ax3.set_ylabel("Frequency",fontsize=18,color=my_grey)
ax3.spines['top'].set_color('none')
ax3.spines['right'].set_color('none')
ax3.spines['bottom'].set_color(my_grey)
ax3.spines['left'].set_color(my_grey)
ax3.tick_params(labelcolor=my_grey,color=my_grey)
leg = ax3.legend(loc="upper right", fontsize=16)
for text in leg.get_texts():
	plt.setp(text, color=my_grey)
plt.savefig(sys.argv[2]+"_dp_qual_mq",dpi=300)
plt.close()


#cSNP hSNP correlation
slope, intercept, r_value, p_value, std_err = stats.linregress(num_csnps,num_hsnps)
print("R^2 between cSNPs and hSNPs",r_value**2)
print('')

# dp info
print('DP INFO')
print('minimum csnp dp value',np.min(csnp_dp))
print('mean read depth amongst cSNPs',np.mean(csnp_dp))
print('standard deviation of cSNP dp values',np.std(csnp_dp))
print('maximum csnp dp value',np.max(csnp_dp))
print('minimum hsnp dp value',np.min(hsnp_dp))
print('mean read depth amongst hSNPs',np.mean(hsnp_dp))
print('standard deviation of hSNP dp values',np.std(hsnp_dp))
print('maximum hsnp dp value',np.max(hsnp_dp))
print('')

# qual info
print('QUAL INFO')
print('minimum base quality amongst csnps',np.min(csnp_qual))
print('mean base quality amongst csnps',np.mean(csnp_qual))
print('standard deviation of cSNP qual values', np.std(csnp_qual))
print('maximum base quality amongst csnps',np.max(csnp_qual))
print('minimum base quality amongst hsnps',np.min(hsnp_qual))
print('mean base quality amongst hsnps',np.mean(hsnp_qual))
print('standard deviation of hSNP qual values',np.std(hsnp_qual))
print('maximum base quality amongst hsnps',np.max(hsnp_qual))
print('')

# mq info
print('MQ INFO')
print('minimum mapping quality amongst csnps',np.min(csnp_qual))
print('mean mapping quality amongst csnps',np.mean(csnp_mq))
print('standard deviation of cSNP mq values',np.std(csnp_mq))
print('maximum mapping quality amongst csnps',np.max(csnp_mq))
print('minimum mapping quality amongst hsnps',np.max(hsnp_mq))
print('mean mapping quality amongst hsnps',np.mean(hsnp_mq))
print('standard deviation of hSNP mq values',np.std(hsnp_mq))
print('maximum mapping quality amongst hsnps',np.max(hsnp_mq))
print('')

# sp info
print('SP INFO')
print('minimum sp score amongst csnps',np.min(csnp_sp))
print('mean sp score amongst csnps',np.mean(csnp_sp))
print('standard deviation of csnp sp values',np.std(csnp_sp))
print('mean hsnp sp score',np.mean(hsnp_sp))
print('std hsnp sp score',np.std(hsnp_sp))
print('mann whitney u sp p-value',stats.mannwhitneyu(csnp_sp,hsnp_sp))
print('95th percentile of csnp sp scores',np.percentile(csnp_sp,95))
print('maximum observed csnp sp score', max(csnp_sp))
print('percentage of hsnp sp scores >95% percentile of csnp_sp scores',   float(sum(i > max(csnp_sp) for i in hsnp_sp))/float(len(hsnp_sp)))
print('maximum observed hsnp sp score',max(hsnp_sp))
print('minimum observed csnp sp score',min(csnp_sp))
print('')

# BQB info
print('BQB INFO')
print('minimum bqb amongst csnps',np.min(csnp_bqb))
print('mean BQB amongst csnps',np.mean(csnp_bqb))
print('standard deviation of csnp bqb values',np.std(csnp_bqb))
print('maximum bqb amongst csnp',np.max(csnp_bqb))
print('minimum bqb amongst hsnps',np.min(hsnp_bqb))
print('mean bqb amongst hsnps',np.mean(hsnp_bqb))
print('standard deviation of hsnp bqb values',np.std(hsnp_bqb))
print('maximum bqb amongst hsnps',np.max(hsnp_bqb))
print('')

#MQB info
print('MBQ INFO')
print('minimum mqb amongst csnps',np.min(csnp_mqb))
print('mean mqb amongst csnps',np.mean(csnp_mqb))
print('standard deviation of csnp mqb values',np.std(csnp_mqb))
print('maximum mqb amongst csnps',np.max(csnp_mqb))
print('minimum mqb amongst hsnps',np.min(hsnp_mqb))
print('mean mqb amongst hsnps',np.mean(hsnp_mqb))
print('standard deviation of hsnp mqb values',np.mean(hsnp_mqb))
print('maximum mqb amongst hsnsp',np.max(hsnp_mqb))
print('')

# rpb info
print('RPB INFO')
print('minimum rpb amongst csnps',np.min(csnp_rpb))
print('mean rpb amongst csnps', np.mean(csnp_rpb))
print('std rpb amongst csnps',np.std(csnp_rpb))
print('maximum rpb amongst csnps',np.max(csnp_rpb))
print('minimum rpb amongst hsnps',np.min(hsnp_rpb))
print('mean rpb amongst hsnps',np.mean(hsnp_rpb))
print('std rpb amongst hsnps',np.std(csnp_rpb))
print('maximum rpb amongst hsnps',np.max(hsnp_rpb))
print('')


# thesis cSNP distance info
print('maximum distance from reference genome',max(num_csnps))
print('minimum pairwise cSNP distance',min(csnp_dist))
print('maximum pairwise cSNP distance',max(csnp_dist))
print('mean pairwise cSNP distance',np.mean(csnp_dist))
print('std pairwise csnp distance',np.std(csnp_dist))

print datetime.now() - startTime
