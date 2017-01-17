#	Mozes Blom, 2 June 2016
#	mozes.blom@gmail.com
# 	The aim of the current script is to generate a mitogenome for each individual from sequence capture data, using a distantly related mitogenome as a reference.
# 	Simple mapping was much less effective and generating an initial reference genome using iterative bating, resulted in  much more effective mapping --> improved mitogenome.
#	One note to still consider, the exact requirements for calling the consensus sequence are not clear. I cross verified for a few samples and was generally satisfied; but an
#	additional analysis that double checks the genotype calling --> read depth would be advisable.

#########################
# Import modules
#########################
import os
import subprocess
import csv
from Bio import SeqIO
from Bio.Seq import *
from natsort import natsorted
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import *
#########################
# Input
#########################

# Main Directory
home_path = "/home/moritz/Desktop/moos/macropus/"

# Specify directory with cleaned paired-end sequence reads (i.e. forward (1), reverse (2) and unpaired (u))
clean_reads_path = "/home/moritz/Desktop/moos/macropus/reads/"

# Specify path to initial species reference, MAKE SURE THAT REFERENCE IS .fa RATHER THAN .fasta
initial_reference_path = "/home/moritz/Desktop/moos/macropus/ref/wallaby_consensus.fa"
ref_strain = "wallaby_ref"

# Specify the list with all individuals to be included in the analysis, names should exactly match the library ID's - read file names (i.e. lib_id_1_final.fastq.gz)
individuals = os.path.join(home_path, '161124_sp_libs.txt')

# Name for output directory, will be created in the 'home' directory specified above
results_dir = "161124_mitogenome"

# Specify the path to each of the executables
interleave_fastq = "/home/moritz/Desktop/moos/bin/MITObim-master/misc_scripts/interleave-fastqgz-MITOBIM.py"
mira =	"/home/moritz/Desktop/moos/bin/mira_4.0.2_linux-gnu_x86_64_static/bin/mira"							# v.4.0.2
mira_dir = "/home/moritz/Desktop/moos/bin/mira_4.0.2_linux-gnu_x86_64_static/bin/"								# For some reason mitobim wants the directory with mira, not the path to the executable
mitobim = "/home/moritz/Desktop/moos/bin/MITObim-master/MITObim_1.8.pl"											# v.1.8
bowtie2_index =	"/home/moritz/Desktop/moos/bin/bowtie2-2.2.9/bowtie2-build"														# v.2.1.0
bowtie2 = "/home/moritz/Desktop/moos/bin/bowtie2-2.2.9/bowtie2"																	# v.2.1.0
samtools = "/home/moritz/Desktop/moos/bin/samtools-1.3.1/samtools"												# v.1.3.1
bcftools = "/home/moritz/Desktop/moos/bin/bcftools-1.3.1/bcftools"												# v.1.3.1
tabix = "/home/moritz/Desktop/moos/bin/htslib-1.3.1/tabix"												# v.1.3.1
# include an aligner!!!


#########################
# Which Analysis to run? 1 = yes, 0 = no
#########################

complete = 1        # Run complete analysis

unzip_reads = 1 			# Unzip the fastq.gz files
indiv_mitobim = 1   		# Run mitobim to generate a mitochondrial reference genome for each individual
indiv_bowtie = 1    		# Map reads to a reference mitogenome (i.e. the mitobim reference) and generate a consensus sequence
collect_consensus = 1     	# Collect all consensus sequences in one fasta file? Here we can also adjust each consensus by coverage if desired
#align = 0 					# Generate a multiple sequence alignment

#########################
# Parameter settings
#########################

# Do the readnames for each library still need a /1 or /2 appended?
add_rp = 0 		# 1 = add a /1 or /2, 0 do nothing
read_id = '@D9S4KXP1'
rm_rp = 0 		# 1 = remove /1 or /2, 0 do nothing

# Mira settings
mira_threads = 2 	# Max. number of hard threads that Mira can use
mira_mem = 8	# Max. amount of RAM that Mira can use

## Use a minimum coverage cut-off? 1 = yes, 0 = no. If so, specify minimum read depth
adjust_coverage = 1
min_read_dp = 3

## Keep FASTQ files (both original and interleaved)? 1 = yes, 0 = no. If yes (1) then the files will be removed
keep_fastq = 0
## Keep MIRA assemblies? 1 = yes, 0 = no. If yes (1) then the folder will be compressed
keep_mira = 0
## Keep MITOBIM assemblies? 1 = yes, 0 = no. If yes (1) then the folder will be compressed
keep_mitobim = 0
## Keep SAM file ? 1 = yes, 0 = no. (strongly advise to remove, very large file!!)
keep_sam = 0
## Keep BAM file ? 1 = yes, 0 = no. (also large file, but potentially useful)
keep_bam = 0
## Keep VCF file ? 1 = yes, 0 = no. (advise to keep)
keep_vcf = 1

#########################
# Functions
#########################

def generate_conf(output_dir, individual, ref_path, ref_strain, reads_path_1_fq, reads_path_2_fq, reads_path_u_fq, threads, max_mem):
	output_path = os.path.join(output_dir, "manifest.conf")
	pf_file = open(output_path, 'w')
   	pf_file.write("project = %s\n" % (individual))
	pf_file.write("job = genome,mapping,accurate\n" )
	pf_file.write("parameters = -GE:not=%s -GE:mps=%s -NW:cnfs=warn -NW:mrnl=0 -NW:cac=no -AS:nop=1 SOLEXA_SETTINGS -CO:msr=no\n" % (str(threads), str(max_mem)))
	pf_file.write("\n")
	pf_file.write("readgroup\n")
	pf_file.write("is_reference\n")
	pf_file.write("data = %s\n" % (ref_path))
	pf_file.write("strain = %s\n" % (ref_strain))
	pf_file.write("\n")
	pf_file.write("readgroup = DataIlluminaPairedLib\n")
	pf_file.write("autopairing\n")
	pf_file.write("data = %s %s\n" % (reads_path_1_fq, reads_path_2_fq))
	pf_file.write("technology = solexa\n")
	pf_file.write("template_size = 50 1000 autorefine\n")
	pf_file.write("segment_placement = ---> <---\n")
	pf_file.write("strain = %s\n" % (individual))
	pf_file.write("\n")
	pf_file.write("readgroup = DataIlluminaUnpaired\n")
	pf_file.write("data = %s\n" % (reads_path_u_fq))
	pf_file.write("technology = solexa\n")
	pf_file.write("strain = %s\n" % (individual))
	pf_file.close()

def add_12_readpair(fq_path, fr, read_id):
	new_fq = []
	with open(fq_path) as fastq:
		for read in fastq:
			if read.startswith(read_id):
				read_string = read.rsplit('\n')[0]
				tmp_read = read_string + "/" + str(fr)
			else:
				tmp_read = read.rsplit('\n')[0]
			new_fq.append(tmp_read)
	return new_fq

def rm_12_readpair(fq_path):
	new_fq = []
	with open(fq_path) as fastq:
		for read in fastq:
			if read.endswith("/1\n"):
				read_string = read.rsplit('/1')[0]
			elif read.endswith("/2\n"):
				read_string = read.rsplit('/2')[0]
			else:
				read_string = read.rsplit('\n')[0]
			new_fq.append(read_string)
	return new_fq

def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		pass
	return False

def run_mira(mira_path, conf_path, work_dir):
	subprocess.call("'%s' '%s' -c '%s'" % (mira_path, conf_path, work_dir), shell=True)

def run_bowtie2_index(bt_index_path, ref_path, index_base):
	subprocess.call("'%s' '%s' '%s'" % (bt_index_path, ref_path, index_base), shell=True)

def run_bowtie2(bt_path, index_base, reads_1_path, reads_2_path, reads_u_path, sam_out_path):
	subprocess.call("'%s' -x '%s' -1 '%s' -2 '%s' -U '%s' -S '%s'" % (bt_path, index_base, reads_1_path, reads_2_path, reads_u_path, sam_out_path), shell=True)

def run_mitobim(mitobim, mira_indiv_strain, mira_ref_strain, reads, maf_path, mira_dir):
	subprocess.call("'%s' -start 1 -end 100 -sample '%s' -ref '%s' -readpool '%s' -maf '%s' --pair --NFS_warn_only --clean --mirapath '%s'" % (mitobim, mira_indiv_strain, mira_ref_strain, reads, maf_path, mira_dir), shell=True)	

def unzip_fq_gz(gz_path, fq_path):
	subprocess.call('gunzip -c %s > %s' % (gz_path, fq_path), shell=True)

def create_dir(dir_path):
	subprocess.call("mkdir '%s'" % (dir_path), shell=True)

def rm(rm_path, file_or_dir):
	# 0 = file
	if file_or_dir == 0:
		subprocess.call("rm '%s'" % (rm_path), shell=True)
	# thus must be folder
	else:
		subprocess.call("rm -r '%s'" % (rm_path), shell=True)

#########################
# Automated Analysis
#########################

# Create an output directory in case it doesn't exist yet
results_dir_path = os.path.join(home_path, results_dir)
if (complete == 1) or (indiv_mitobim == 1):
    if os.path.isdir(results_dir_path):
        sys.exit("The results folder: %s ALREADY EXISTS, please specify a new folder name" % (results_dir_path))
    else:
        create_dir(results_dir_path)
else:
    pass

# Create a list of all target individuals
indivs_list = []
with open(individuals) as id_list:
	for line in id_list:
		name = line.rsplit('\n')[0]
		indivs_list.append(name)

# unzip reads (add /1 or /2 if needed)
if (complete == 1) or (unzip_reads == 1):
	for name in indivs_list:
		for lib in ["1", "2", "u"]:
			gz_path = os.path.join(clean_reads_path, (name + "_" + lib + "_final.fastq.gz"))
			fq_path = os.path.join(clean_reads_path, (name + "_" + lib + "_final.fastq"))
			unzip_fq_gz(gz_path, fq_path)
			if add_rp == 1:
			# Comment the following section out in case the /1 or /2 are already appended to the read id in the fastq files
				if (lib == "1") or (lib == "2"):
					fq_updated = add_12_readpair(fq_path, lib, read_id)
					fq_final = []
					record_counter = 0
					while record_counter < len(fq_updated):
						record = fq_updated[record_counter]
						if record.startswith(read_id):
							if (is_number(record[-3]) == True) and (is_number(record[-4]) == True) and (is_number(record[-5]) == True):
								fq_final.append(record + '\n')
								record_counter += 1
							else:
								print record
								record_counter += 4
						else:
							fq_final.append(record + '\n')
							record_counter += 1
					with open(fq_path, "w") as fq_file:
						fq_file.writelines(fq_final)
				else:
					with open(fq_path) as fastq:
						fq_updated = []
						for read in fastq:
							fq_updated.append(read.rsplit('\n')[0])
					fq_final = []
					record_counter = 0
					while record_counter < len(fq_updated):
						record = fq_updated[record_counter]
						if record.startswith(read_id):
							if (is_number(record[-1]) == True) and (is_number(record[-2]) == True) and (is_number(record[-3]) == True):
								fq_final.append(record + '\n')
								record_counter += 1
							else:
								print record
								record_counter += 4
						else:
							fq_final.append(record + '\n')
							record_counter += 1
					with open(fq_path, "w") as fq_file:
						fq_file.writelines(fq_final)
			elif rm_rp == 1:
				fq_updated = rm_12_readpair(fq_path)
				fq_final = []
				record_counter = 0
				while record_counter < len(fq_updated):
					record = fq_updated[record_counter]
					if record.startswith(read_id):
						if (is_number(record[-1]) == True) and (is_number(record[-2]) == True) and (is_number(record[-3]) == True):
							fq_final.append(record + '\n')
							record_counter += 1
						else:
							print record
							record_counter += 4
					else:
						fq_final.append(record + '\n')
						record_counter += 1
				with open(fq_path, "w") as fq_file:
					fq_file.writelines(fq_final)
			# ASSUMING THE /1 AND /2 ARE ALREADY APPENDED TO EACH READ NAME!!!!
			else:	
				with open(fq_path) as fastq:
					fq_updated = []
					for read in fastq:
						fq_updated.append(read.rsplit('\n')[0])
				fq_final = []
				record_counter = 0
				while record_counter < len(fq_updated):
					record = fq_updated[record_counter]
					if record.startswith(read_id):
						if (is_number(record[-3]) == True) and (is_number(record[-4]) == True) and (is_number(record[-5]) == True):
							fq_final.append(record + '\n')
							record_counter += 1
						else:
							print record
							record_counter += 4
					else:
						fq_final.append(record + '\n')
						record_counter += 1
				with open(fq_path, "w") as fq_file:
					fq_file.writelines(fq_final)
else:
	pass

mb_path = os.path.join(results_dir_path, '1.mitobim')
# Generate reference sequences for each individual using mitoBIM
if (complete == 1) or (indiv_mitobim == 1):
	create_dir(mb_path)
	for indiv in indivs_list:
		create_dir(os.path.join(mb_path, indiv))
		# generate config file for mira: generate_conf(output folder, individual, reference fasta, ref. fasta strain, clean reads 1/2/u)
		generate_conf(os.path.join(mb_path, indiv), indiv, initial_reference_path, ref_strain, os.path.join(clean_reads_path, (indiv + "_1_final.fastq")), os.path.join(clean_reads_path, (indiv + "_2_final.fastq")), os.path.join(clean_reads_path, (indiv + "_u_final.fastq")), mira_threads, mira_mem)
		conf_file_path = os.path.join(mb_path, indiv, 'manifest.conf')
		# run mira (mira path, config file, output directory)
		run_mira(mira, conf_file_path, os.path.join(mb_path, indiv))
		# generate an interleaved fastq file for the /1 and /2 readpairs
		fq_interleaved_path = os.path.join(mb_path, indiv, (indiv + "_interleaved.fastq"))
		subprocess.call('%s %s %s > %s' % (interleave_fastq, os.path.join(clean_reads_path, (indiv + "_1_final.fastq")), os.path.join(clean_reads_path, (indiv + "_2_final.fastq")), fq_interleaved_path), shell=True)
		# run MITObim
		maf = os.path.join(mb_path, indiv, (indiv + '_assembly'), (indiv + '_d_results'), (indiv + '_out.maf'))
		os.chdir(os.path.join(mb_path, indiv))
		run_mitobim(mitobim, indiv, ref_strain, fq_interleaved_path, maf, mira_dir)
		if keep_fastq == 0:
			rm(os.path.join(mb_path, (indiv + "_interleaved.fastq")), 0)
		else:
			pass
		if keep_mira == 0:
			rm(os.path.join(mb_path, indiv, (indiv + '_assembly')), 1)
		else:
			pass


bt_path = os.path.join(results_dir_path, '2.bowtie2')
# Map reads for each individual back to it's own MITObim reference
if (complete == 1) or (indiv_bowtie == 1):
	create_dir(bt_path)
	for indiv in indivs_list:
		create_dir(os.path.join(bt_path, indiv))
		# Need to find out which was the final MITObim iteration
		iter_list = []
		# This Try statement was included, because the os.listdir call did not work if there was a .XX file in the folder.
		try:
			subprocess.call("rm %s" % (os.path.join(mb_path, indiv, ".DS_store")))
		except:
			pass
		try:
			for mb_iter in os.listdir(os.path.join(mb_path, indiv)):
				if mb_iter.startswith('iteration'):
					iter_list.append(mb_iter)
				else:
					pass
			iter_list = natsorted(iter_list)
			final_iter = iter_list[-1]
			# Copy MITObim reference to Bowtie folder
			subprocess.call("cp %s %s" % (os.path.join(mb_path, indiv, final_iter, (indiv + "-" + ref_strain + "_assembly"), (indiv + "-" + ref_strain + "_d_results"), (indiv + "-" + ref_strain + "_out_" + ref_strain + ".unpadded.fasta")), os.path.join(bt_path, indiv, (indiv + "_mitobim_ref.fasta"))), shell=True)
			os.chdir(os.path.join(bt_path, indiv))
			# Run Bowtie2
			run_bowtie2_index(bowtie2_index, os.path.join(bt_path, indiv, (indiv + "_mitobim_ref.fasta")), (indiv + "_ref"))
			run_bowtie2(bowtie2, (indiv + "_ref"), os.path.join(clean_reads_path, (indiv + "_1_final.fastq")), os.path.join(clean_reads_path, (indiv + "_2_final.fastq")), os.path.join(clean_reads_path, (indiv + "_u_final.fastq")), os.path.join(bt_path, indiv, (indiv + "_mapped_reads.sam")))
			# Convert the SAM file into a (sorted) BAM file
			subprocess.call("%s view -bS %s > %s" % (samtools, os.path.join(bt_path, indiv, (indiv + "_mapped_reads.sam")), os.path.join(bt_path, indiv, (indiv + "_mapped_reads.bam"))), shell=True)
			# Remove unsorted BAM file to free up space
			if keep_sam == 0:
				rm(os.path.join(bt_path, indiv, (indiv + "_mapped_reads.sam")), 0)
			else:
				pass
			subprocess.call("%s sort %s > %s" % (samtools, os.path.join(bt_path, indiv, (indiv + "_mapped_reads.bam")), os.path.join(bt_path, indiv, (indiv + "_mapped_reads_sorted.bam"))), shell=True)
			# Remove unsorted BAM file to free up space
			if keep_bam == 0:
				rm(os.path.join(bt_path, indiv, (indiv + "_mapped_reads.bam")), 0)
			else:
				pass
			# Calculate read depth
			subprocess.call("%s depth -a -q 10 -Q 10 %s > %s" % (samtools, os.path.join(bt_path, indiv, (indiv + "_mapped_reads_sorted.bam")), os.path.join(bt_path, indiv, (indiv + "_mapped_reads_depth.txt"))), shell=True)
			# Calculate genotype likelihoods
			subprocess.call("%s mpileup -Ou -f %s %s| %s call -mv -Oz -o %s" % (samtools, os.path.join(bt_path, indiv, (indiv + "_mitobim_ref.fasta")), os.path.join(bt_path, indiv, (indiv + "_mapped_reads_sorted.bam")), bcftools, os.path.join(bt_path, indiv, (indiv + "_calls.vcf.gz"))), shell=True)
			subprocess.call("%s %s" % (tabix, os.path.join(bt_path, indiv, (indiv + "_calls.vcf.gz"))), shell=True)
			# Generate a consensus sequence, incorporating the genotype calls from the vcf file
			subprocess.call("cat %s | %s consensus %s > %s" % (os.path.join(bt_path, indiv, (indiv + "_mitobim_ref.fasta")), bcftools, os.path.join(bt_path, indiv, (indiv + "_calls.vcf.gz")), os.path.join(bt_path, indiv, (indiv + "_consensus.fa"))), shell=True)
		except:
			print ('DID NOT MAP READS TO INDIVIDUAL %s MOST LIKELY SINCE NO MITOBIM ITERATIONS WERE FOUND' % indiv)


cs_path = os.path.join(results_dir_path, '3.mitogenomes')
# Gather all inferred consensus sequences, for each individual, and adjust by coverage if desired
if (complete == 1) or (collect_consensus == 1):
	create_dir(cs_path)
	all_consens = []
	for indiv in indivs_list:
		create_dir(os.path.join(cs_path, indiv))
		os.chdir(os.path.join(cs_path, indiv))
		try:
			subprocess.call("cp %s %s" % (os.path.join(bt_path, indiv, (indiv + "_consensus.fa")), os.path.join(cs_path, indiv)), shell=True)
			if adjust_coverage == 1:
				run_bowtie2_index(bowtie2_index, os.path.join(cs_path, indiv, (indiv + "_consensus.fa")), (indiv + "_ref"))
				run_bowtie2(bowtie2, (indiv + "_ref"), os.path.join(clean_reads_path, (indiv + "_1_final.fastq")), os.path.join(clean_reads_path, (indiv + "_2_final.fastq")), os.path.join(clean_reads_path, (indiv + "_u_final.fastq")), os.path.join(cs_path, indiv, (indiv + "_mapped_reads.sam")))
				# Remove fastq files to free up space
				if keep_fastq == 0:
					for lib in ["1", "2", "u"]:
						rm(os.path.join(clean_reads_path, (indiv + "_" + lib + "_final.fastq")), 0)
				else:
					pass
				# Convert the SAM file into a (sorted) BAM file
				subprocess.call("%s view -bS %s > %s" % (samtools, os.path.join(cs_path, indiv, (indiv + "_mapped_reads.sam")), os.path.join(cs_path, indiv, (indiv + "_mapped_reads.bam"))), shell=True)
				# Remove unsorted BAM file to free up space
				if keep_sam == 0:
					rm(os.path.join(cs_path, indiv, (indiv + "_mapped_reads.sam")), 0)
				else:
					pass
				subprocess.call("%s sort %s > %s" % (samtools, os.path.join(cs_path, indiv, (indiv + "_mapped_reads.bam")), os.path.join(cs_path, indiv, (indiv + "_mapped_reads_sorted.bam"))), shell=True)
				# Remove unsorted BAM file to free up space
				if keep_bam == 0:
					rm(os.path.join(cs_path, indiv, (indiv + "_mapped_reads.bam")), 0)
				else:
					pass
				# Calculate read depth
				subprocess.call("%s depth -a -q 10 -Q 10 %s > %s" % (samtools, os.path.join(cs_path, indiv, (indiv + "_mapped_reads_sorted.bam")), os.path.join(cs_path, indiv, (indiv + "_mapped_reads_depth.txt"))), shell=True)
				indiv_cons_fa = SeqIO.read(os.path.join(cs_path, indiv, (indiv + "_consensus.fa")), 'fasta')
				raw_cons = indiv_cons_fa.seq
				adj_cons = ''
				with open(os.path.join(cs_path, indiv, (indiv + "_mapped_reads_depth.txt"))) as cov_stats:
					depth_stats = csv.DictReader(cov_stats, delimiter = '\t', fieldnames = ['ref', 'pos', 'depth'])
					for row in depth_stats:
						if int(row['depth']) >= min_read_dp:
							adj_cons = adj_cons + raw_cons[int(row['pos'])]
						else:
							adj_cons = adj_cons + 'n'
				record = SeqRecord(Seq(adj_cons, generic_alphabet), id=indiv, description=indiv)
				all_consens.append(record)
				output_handle = open(os.path.join(cs_path, indiv, (indiv + "_consensus.fa")), "w")
				SeqIO.write(record, output_handle, "fasta")
			else:
				indiv_cons_fa = SeqIO.read(os.path.join(cs_path, indiv, (indiv + "_consensus.fa")), 'fasta')
				record = SeqRecord(indiv_cons_fa.seq, id=indiv, description=indiv)
				all_consens.append(record)
				output_handle = open(os.path.join(cs_path, indiv, (indiv + "_consensus.fa")), "w")
				SeqIO.write(record, output_handle, "fasta")
		except:
			pass
	output_handle = open(os.path.join(cs_path, ("all_consensus.fa")), "w")
	SeqIO.write(all_consens, output_handle, "fasta")
	output_handle.close()










