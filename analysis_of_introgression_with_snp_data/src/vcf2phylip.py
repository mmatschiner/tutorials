#!/usr/bin/env python2
# -*- coding: utf-8 -*-



'''
The script converts a collection of SNPs in VCF format into a PHYLIP, FASTA, 
NEXUS, or binary NEXUS file for phylogenetic analysis. The code is optimized
to process VCF files with sizes >1GB. For small VCF files the algorithm slows
down as the number of taxa increases (but is still fast).
Any ploidy is allowed, but binary NEXUS is produced only for diploid VCFs.
'''



__author__      = "Edgardo M. Ortiz"
__credits__     = "Juan D. Palacio-Mejía"
__version__     = "2.0"
__email__       = "e.ortiz.v@gmail.com"
__date__        = "2019-01-15"



import sys
import os
import argparse
import gzip



def main():
	parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument("-i", "--input", action="store", dest="filename", required=True,
		help="Name of the input VCF file, can be gzipped")
	parser.add_argument("-m", "--min-samples-locus", action="store", dest="min_samples_locus", type=int, default=4,
		help="Minimum of samples required to be present at a locus, default=4 since is the minimum for phylogenetics.")
	parser.add_argument("-o", "--outgroup", action="store", dest="outgroup", default="",
		help="Name of the outgroup in the matrix. Sequence will be written as first taxon in the alignment.")
	parser.add_argument("-p", "--phylip-disable", action="store_true", dest="phylipdisable", default=False,
		help="A PHYLIP matrix is written by default unless you enable this flag")
	parser.add_argument("-f", "--fasta", action="store_true", dest="fasta", default=False,
		help="Write a FASTA matrix, disabled by default")
	parser.add_argument("-n", "--nexus", action="store_true", dest="nexus", default=False,
		help="Write a NEXUS matrix, disabled by default")
	parser.add_argument("-b", "--nexus-binary", action="store_true", dest="nexusbin", default=False,
		help="Write a binary NEXUS matrix for analysis of biallelic SNPs in SNAPP, disabled by default")
	args = parser.parse_args()



	filename = args.filename

    # Prepare the opener if the SAM file is compressed
	if filename.endswith(".gz"):
		opener = gzip.open
	else:
		opener = open

	min_samples_locus = args.min_samples_locus
	outgroup = args.outgroup
	phylipdisable = args.phylipdisable
	fasta = args.fasta
	nexus = args.nexus
	nexusbin = args.nexusbin



	# Dictionary of IUPAC ambiguities for nucleotides
	# '*' means deletion for GATK (and other software?)
	# Deletions are ignored when making the consensus
	# Dictionary to translate IUPAC ambiguities, lowercase letters are used when "*" or "N" were present for a position,
	# however, software like Genious for example are case insensitive and will imply ignore capitalization
	amb = {"*":"-", "A":"A", "C":"C", "G":"G", "N":"N", "T":"T",
		   "*A":"a", "*C":"c", "*G":"g", "*N":"n", "*T":"t",
		   "AC":"M", "AG":"R", "AN":"a", "AT":"W", "CG":"S",
		   "CN":"c", "CT":"Y", "GN":"g", "GT":"K", "NT":"t",
		   "*AC":"m", "*AG":"r", "*AN":"a", "*AT":"w", "*CG":"s",
		   "*CN":"c", "*CT":"y", "*GN":"g", "*GT":"k", "*NT":"t",
		   "ACG":"V", "ACN":"m", "ACT":"H", "AGN":"r", "AGT":"D",
		   "ANT":"w", "CGN":"s", "CGT":"B", "CNT":"y", "GNT":"k",
		   "*ACG":"v", "*ACN":"m", "*ACT":"h", "*AGN":"r", "*AGT":"d",
		   "*ANT":"w", "*CGN":"s", "*CGT":"b", "*CNT":"y", "*GNT":"k",
		   "ACGN":"v", "ACGT":"N", "ACNT":"h", "AGNT":"d", "CGNT":"b",
		   "*ACGN":"v", "*ACGT":"N", "*ACNT":"h", "*AGNT":"d", "*CGNT":"b",
		   "*ACGNT":"N"}



	# Dictionary for translating biallelic SNPs into SNAPP, only for diploid VCF
	# 0 is homozygous reference
	# 1 is heterozygous
	# 2 is homozygous alternative
	gen_bin = {"./.":"?",
			   ".|.":"?",
			   "0/0":"0",
			   "0|0":"0",
			   "0/1":"1",
			   "0|1":"1",
			   "1/0":"1",
			   "1|0":"1",
			   "1/1":"2",
			   "1|1":"2"}



	############################
	# Process header of VCF file
	ploidy = 0
	gt_idx = []
	missing = ""

	with opener(filename) as vcf:

		# Create a list to store sample names
		sample_names = []

		# Keep track of longest sequence name for padding with spaces in the output file
		len_longest_name = 0

		# Look for the line in the VCF header with the sample names
		for line in vcf:
			if line.startswith("#CHROM"):

				# Split line into fields
				broken = line.strip("\n").split("\t")

				# If the minimum-samples-per-locus parameter is larger than the number of
				# species in the alignment make it the same as the number of species
				if min_samples_locus > len(broken[9:]):
					min_samples_locus = len(broken[9:])

				# Create a list of sample names and the keep track of the longest name length
				for i in range(9, len(broken)):
					name_sample = broken[i].replace("./","") # GATK adds "./" to sample names sometimes
					sample_names.append(name_sample)
					len_longest_name = max(len_longest_name, len(name_sample))

			# Find out the ploidy of the genotypes, just distinguishes if sample is not haploid vs n-ploid
			elif not line.startswith("#") and line.strip("\n") != "":
				while ploidy == 0 and missing == "":
					broken = line.strip("\n").split("\t")
					for j in range(9, len(broken)):
						if ploidy == 0:
							if broken[j].split(":")[0][0] != ".":
								ploidy = (len(broken[j].split(":")[0]) // 2) + 1
								gt_idx = [i for i in range(0, len(broken[j].split(":")[0]), 2)]
								missing = "/".join(["." for i in range(len(gt_idx))])
								# print missing
								print len(broken)
								print broken[j]
								print gt_idx
								# print ploidy
								break
				
	vcf.close()

	print("\nConverting file " + filename + ":\n")
	print("Number of samples in VCF: " + str(len(sample_names)))


	####################
	# SETUP OUTPUT FILES 
	# Output filename will be the same as input file, indicating the minimum of samples specified
	if filename.endswith(".gz"):
		outfile = filename.replace(".vcf.gz",".min"+str(min_samples_locus))
	else:
		outfile = filename.replace(".vcf",".min"+str(min_samples_locus))

	# We need to create an intermediate file to hold the sequence data 
	# vertically and then transpose it to create the matrices
	if fasta or nexus or not phylipdisable:
		temporal = open(outfile+".tmp", "w")
	
	# if binary NEXUS is selected also create a separate temporal
	if nexusbin:
		if ploidy == 2:
			temporalbin = open(outfile+".bin.tmp", "w")
		else:
			print("Binary NEXUS not available for "+str(ploidy)+"-ploid VCF.")
			nexusbin = False



	##################
	# PROCESS VCF FILE
	index_last_sample = len(sample_names)+9

	# Start processing SNPs of VCF file
	with opener(filename) as vcf:

		# Initialize line counter
		snp_num = 0
		snp_accepted = 0
		snp_shallow = 0
		snp_multinuc = 0
		snp_biallelic = 0
		while 1:

			# Load large chunks of file into memory
			vcf_chunk = vcf.readlines(50000)
			if not vcf_chunk:
				break

			# Now process the SNPs one by one
			for line in vcf_chunk:
				if not line.startswith("#") and line.strip("\n") != "": # pyrad sometimes produces an empty line after the #CHROM line

					# Split line into columns
					broken = line.strip("\n").split("\t")
					for g in range(9,len(broken)):
						if broken[g].split(":")[0][0] == ".":
							broken[g] = missing

					# Keep track of number of genotypes processed
					snp_num += 1

					# Print progress every 500000 lines
					if snp_num % 500000 == 0:
						print(str(snp_num)+" genotypes processed.")

					# Check if the SNP has the minimum of samples required
					if (len(broken[9:]) - ''.join(broken[9:]).count(missing)) >= min_samples_locus:
						
						# Check that ref genotype is a single nucleotide and alternative genotypes are single nucleotides
						if len(broken[3]) == 1 and (len(broken[4])-broken[4].count(",")) == (broken[4].count(",")+1):

							# Add to running sum of accepted SNPs
							snp_accepted += 1

							# If nucleotide matrices are requested
							if fasta or nexus or not phylipdisable:

								# Create a dictionary for genotype to nucleotide translation
								# each SNP may code the nucleotides in a different manner
								nuc = {str(0):broken[3], ".":"N"}
								for n in range(len(broken[4].split(","))):
									nuc[str(n+1)] = broken[4].split(",")[n]

								# Translate genotypes into nucleotides and the obtain the IUPAC ambiguity
								# for heterozygous SNPs, and append to DNA sequence of each sample
								#print index_last_sample
								#print nuc[broken[9][0]]
                 						site_tmp = ''.join([amb[''.join(sorted(set([nuc[broken[i][j]] for j in gt_idx])))] for i in range(9, index_last_sample)])

								# Write entire row of single nucleotide genotypes to temporary file
								temporal.write(site_tmp+"\n")

								# print nuc
								# print [i.split(":")[0] for i in broken[9:]]
								# print site_tmp

							# Write binary NEXUS for SNAPP if requested
							if nexusbin:

								# Check taht the SNP only has two alleles
								if len(broken[4]) == 1:
									
									# Add to running sum of biallelic SNPs
									snp_biallelic += 1

									# Translate genotype into 0 for homozygous ref, 1 for heterozygous, and 2 for homozygous alt
									binsite_tmp = ''.join([(gen_bin[broken[i][0:3]]) for i in range(9, index_last_sample)])

									# Write entire row to temporary file
									temporalbin.write(binsite_tmp+"\n")

						else:
							# Keep track of loci rejected due to multinucleotide genotypes
							snp_multinuc += 1
							# Keep track of loci rejected due to exceeded missing data
							snp_shallow += 1

					else:
						# Keep track of loci rejected due to exceeded missing data
						snp_shallow += 1

		# Print useful information about filtering of SNPs
		print("Total of genotypes processed: " + str(snp_num))
		print("Genotypes excluded because they exceeded the amount of missing data allowed: " + str(snp_shallow))
		print("Genotypes that passed missing data filter but were excluded for not being SNPs: " + str(snp_multinuc))
		print("SNPs that passed the filters: " + str(snp_accepted)) 
		if nexusbin:
			print("Biallelic SNPs selected for binary NEXUS: " + str(snp_biallelic))
		print("")

	vcf.close()
	if fasta or nexus or not phylipdisable:
		temporal.close()
	if nexusbin:
		temporalbin.close()



	#######################
	# WRITE OUTPUT MATRICES

	if not phylipdisable:
		output_phy = open(outfile+".phy", "w")
		header_phy = str(len(sample_names))+" "+str(snp_accepted)+"\n"
		output_phy.write(header_phy)

	if fasta:
		output_fas = open(outfile+".fasta", "w")

	if nexus:
		output_nex = open(outfile+".nexus", "w")
		header_nex = "#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX=" + str(len(sample_names)) + " NCHAR=" + str(snp_accepted) + ";\n\tFORMAT DATATYPE=DNA" + " MISSING=N" + " GAP=- ;\nMATRIX\n"
		output_nex.write(header_nex)

	if nexusbin:
		output_nexbin = open(outfile+".bin.nexus", "w")
		header_nexbin = "#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX=" + str(len(sample_names)) + " NCHAR=" + str(snp_biallelic) + ";\n\tFORMAT DATATYPE=SNP" + " MISSING=?" + " GAP=- ;\nMATRIX\n"
		output_nexbin.write(header_nexbin)


	# Store index of outgroup in list of sample names
	idx_outgroup = "NA"

	# Write outgroup as first sequence in alignment if the name is specified
	if outgroup in sample_names:
		idx_outgroup = sample_names.index(outgroup)

		if fasta or nexus or not phylipdisable:
			with open(outfile+".tmp") as tmp_seq:
				seqout = ""

				# This is where the transposing happens
				for line in tmp_seq:
					seqout += line[idx_outgroup]

				# Write FASTA line
				if fasta:
					output_fas.write(">"+sample_names[idx_outgroup]+"\n"+seqout+"\n")
				
				# Pad sequences names and write PHYLIP or NEXUS lines
				padding = (len_longest_name + 3 - len(sample_names[idx_outgroup])) * " "
				if not phylipdisable:
					output_phy.write(sample_names[idx_outgroup]+padding+seqout+"\n")
				if nexus:
					output_nex.write(sample_names[idx_outgroup]+padding+seqout+"\n")

				# Print current progress
				print("Outgroup, "+outgroup+", added to the matrix(ces).")

		if nexusbin:
			with open(outfile+".bin.tmp") as bin_tmp_seq:
				seqout = ""

				# This is where the transposing happens
				for line in bin_tmp_seq:
					seqout += line[idx_outgroup]

				# Write line of binary SNPs to NEXUS
				padding = (len_longest_name + 3 - len(sample_names[idx_outgroup])) * " "
				output_nexbin.write(sample_names[idx_outgroup]+padding+seqout+"\n")

				# Print current progress
				print("Outgroup, "+outgroup+", added to the binary matrix.")


	# Write sequences of the ingroup
	for s in range(0, len(sample_names)):
		if s != idx_outgroup:
			if fasta or nexus or not phylipdisable:
				with open(outfile+".tmp") as tmp_seq:
					seqout = ""

					# This is where the transposing happens
					for line in tmp_seq:
						seqout += line[s]

					# Write FASTA line
					if fasta:
						output_fas.write(">"+sample_names[s]+"\n"+seqout+"\n")
					
					# Pad sequences names and write PHYLIP or NEXUS lines
					padding = (len_longest_name + 3 - len(sample_names[s])) * " "
					if not phylipdisable:
						output_phy.write(sample_names[s]+padding+seqout+"\n")
					if nexus:
						output_nex.write(sample_names[s]+padding+seqout+"\n")

					# Print current progress
					print("Sample "+str(s+1)+" of "+str(len(sample_names))+", "+sample_names[s]+", added to the nucleotide matrix(ces).")

			if nexusbin:
				with open(outfile+".bin.tmp") as bin_tmp_seq:
					seqout = ""

					# This is where the transposing happens
					for line in bin_tmp_seq:
						seqout += line[s]

					# Write line of binary SNPs to NEXUS
					padding = (len_longest_name + 3 - len(sample_names[s])) * " "
					output_nexbin.write(sample_names[s]+padding+seqout+"\n")

					# Print current progress
					print("Sample "+str(s+1)+" of "+str(len(sample_names))+", "+sample_names[s]+", added to the binary matrix.")


	if not phylipdisable:
		output_phy.close()
	if fasta:
		output_fas.close()
	if nexus:
		output_nex.write(";\nEND;\n")
		output_nex.close()
	if nexusbin:
		output_nexbin.write(";\nEND;\n")
		output_nexbin.close()

	if fasta or nexus or not phylipdisable:
		os.remove(outfile+".tmp")
	if nexusbin:
		os.remove(outfile+".bin.tmp")


	print( "\nDone!\n")


if __name__ == "__main__":
    main()
