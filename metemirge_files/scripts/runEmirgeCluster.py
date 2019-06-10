#!/usr/bin/env python

"""
Shalabh Thakur
March, 2018
Sycuro Lab, University of Calgary
"""
"""
Description: This script is the main controller from running EMIRGE.
"""
import sys
import os
import glob
import argparse
import re
import shutil
import logging
import datetime
from argparse import RawTextHelpFormatter
from Bio import SeqIO
import statistics


##### MAIN FUNCTION #######
def main():

	"""
	Main Function
	"""
	now = datetime.datetime.now()
	
	c_args = parse_args(__file__)
	
	project_dir=os.path.abspath(c_args['project_dir'])
	config_file=os.path.abspath(c_args['config_file'])
	output_dir=os.path.abspath(c_args['output_dir'])
	data_dir=os.path.abspath(c_args['data_dir'])
	cpu=int(c_args['cpu'])
	sample_name=c_args['sample']
	

	###### Read Config File Function #####
	config_dict={}
	config_parameters={}
	program_name=""
	project_name=""
	
	if(c_args['config_file']):
		config_file=os.path.abspath(c_args['config_file'])
		config_parameters=read_config_file(config_file)
		program_name=config_parameters['PROGRAM']
		project_name=config_parameters['PROJECT_NAME']
		config_dict[program_name]=config_parameters

	###### Get Sample Dir Names ########
	#sample_name=get_sample_name(jobindex,jobindex_file)

	##### get list of sample dir names ####
	sample_dir=os.path.join(data_dir,sample_name)
	
        ##### Define Project Path#####
	project_sample_path=os.path.join(project_dir,sample_name)
        project_path=os.path.join(project_dir,sample_name,project_name)
        
        if not os.path.exists(project_sample_path):
		os.makedirs(project_sample_path)
		

        if not os.path.exists(project_path):
		os.makedirs(project_path)
		
	log_file=os.path.join(project_path,sample_name+"_"+project_name.lower()+"_run.log")
	logging.basicConfig(format="[%(asctime)s] %(levelname)-8s %(message)s",datefmt='%a, %d %b %Y %H:%M:%S',level=logging.INFO)
	logging.info("Logging information for runEmirgeCluster.py run for {}\n".format(sample_name))

	logging.info("Creating sub-directory with sample name {} in main project directory at {}".format(sample_name, project_sample_path))
	logging.info("Creating sub-directory with project name {} for sample {} at {}".format(project_name,sample_name,project_path))

	logging.info("Program Name: {}".format(program_name))
	logging.info("Project Name is {}".format(project_name))
	logging.info("Sample Name is {}".format(sample_name))
	logging.info("Reading sample data directory at {}".format(sample_dir))
	logging.info("Reading config file at {}".format(config_file))
	logging.info("Main Project directory is located at {}".format(project_dir))

	####### Run EMIRGE ####
	if(program_name=="EMIRGE"):
		logging.info("Reading parameters for running EMIRGE program")
		emirge_parameter=config_dict["EMIRGE"]
		logging.info("Read {} parameters from the config file at {}".format(len(emirge_parameter),config_file))

		read_type=emirge_parameter['READ_TYPE']
		read_format=emirge_parameter['READ_FORMAT']
		suffix_read1=emirge_parameter['READ1_SUFFIX']+"."+read_format
		suffix_read2=emirge_parameter['READ2_SUFFIX']+"."+read_format
	
		##### Full path of read files ####
		read_1_file=""
		read_2_file=""

		#### read file names in the sample directories ####
		list_sample_file=get_file_name_from_dir(sample_dir)
		
		logging.info("Reading sequencing read files for sample {} from {}".format(sample_name,sample_dir))
		logging.info("Read type: {}".format(read_type))
		logging.info("Read format: {}\n".format(read_format))

		##### If reads are paired-ended ######
		if read_type=="paired":
			#### If forward and reverse reads are in a seperate file #######
			for read_file in list_sample_file:
				##### If read file has suffix for read1 ####
				if re.search(suffix_read1,read_file):
					read_1_file=os.path.join(sample_dir,read_file)
					logging.info("Suffix used for forward reads: {}".format(suffix_read1))
					logging.info("Reading forward read file {}".format(read_1_file))
					
				elif re.search(suffix_read2,read_file):
				#### If read file has suffix for read 2 ####
					read_2_file=os.path.join(sample_dir,read_file)
					logging.info("Suffix used for reverse reads: {}".format(suffix_read2))
					logging.info("Reading reverse read file {}".format(read_2_file))

		##### If reads are single-ended ######
		if read_type=="single":
			read_1_file=os.path.join(sample_dir,list_sample_file[0])
			logging.info("Reading single read file {}".format(read_1_file))

		print read_1_file
		print read_2_file	
		##### Call function for specific program ######
		logging.info("Running function for executing EMIRGE program on sample {}".format(sample_name))		
		run_emirge(sample_name,read_1_file,read_2_file,project_path,config_dict,cpu)

		logging.info("Running function for renaming EMIRGE-generated fasta file for sample {}".format(sample_name))
		rename_emirge_fasta(sample_name,project_path)

	######## RUN EMIRGE BLAST #####
	elif(program_name=="EMIRGE_BLAST"):

		logging.info("Reading parameters for running nucleotide BLAST program")
		emirge_blast_parameters=config_dict["EMIRGE_BLAST"]
	
		logging.info("Running function for executing nucleotide BLAST program on sample {}".format(sample_name))
		run_emirge_blast(sample_name,sample_dir,project_name,project_path,emirge_blast_parameters)

	elif (program_name=="EMIRGE_SELECT"):

		logging.info("Reading parameters for selecting best blast-hit")
		emirge_select_parameters=config_dict["EMIRGE_SELECT"]

		logging.info("Running function for selecting best blast-hit for sample {}".format(sample_name))
		run_emirge_select(sample_name,sample_dir,project_name,project_path,emirge_select_parameters)

		

##### FUNCTIONS #######
def parse_args(desc):

	"""
	Command line option parser
	
	:param desc: Short program description.
	:type desc: str
	
	:return arg: dict of command line arguments with key=command line
        argument and val=argument value
        
	:return_type: dict
	"""

	parser = argparse.ArgumentParser(description=desc, formatter_class=RawTextHelpFormatter)
	parser._optionals.title="Arguments"
	parser.add_argument("--config_file", help="name and path of the single config file", required=True)
	parser.add_argument("--project_dir", help="name and path of the project", required=True)
	parser.add_argument("--data_dir", help="path to the directory containing samples", required=True)
	parser.add_argument("--output_dir", help="path to the output directory", required=True)
	parser.add_argument("--cpu",type=int, default=4, help="specify number of cpu cores to use\n", required=True)
	parser.add_argument("--sample", help="specify sample name\n", required=True)
	
	args=parser.parse_args()

	return vars(args)

#### Function to read config file ####
def read_config_file(config_file):	

	"""
	This functions reads config files
	"""	
	config_param={}

	with open(config_file,"r") as config:

		for line in config:

			if re.search("#",line):
				continue
			elif re.search('^\s+$',line):
				continue

			line=line.rstrip("\n")
			param_column=line.split(":")
			parameter_name=param_column[0].strip()
			parameter_value=param_column[1].strip()
			print "{}: {}\n".format(parameter_name,parameter_value)
			config_param[parameter_name]=parameter_value		

	return(config_param)


#### Function to read sample dir names #####

def get_dir_list(directory):
	"""
	Input: directory name
	"""
	return [name for name in os.listdir(directory)
		if os.path.isdir(os.path.join(directory, name))]

#### Function to read sample file names ####
def get_file_name_from_dir(directory):
	"""
	Input: directory name
	"""
	return [name for name in os.listdir(directory)
		if os.path.isfile(os.path.join(directory, name))]
	

#### Function to run EMIRGE #####

def run_emirge(sample_name,read_1_file,read_2_file,project_path,config_dict,cpu):
	
	logging.info("Analysis Step: Marker sequence reconstruction")
	emirge_parameter=config_dict["EMIRGE"]
	
	#### create output directory for the sample under processing ####
	sample_out_dir=os.path.join(project_path,"EMIRGE_RECONSTRUCT")

	if not os.path.exists(sample_out_dir):
		os.makedirs(sample_out_dir)
		logging.info("Creating output sub-directory {}".format(sample_out_dir))
	else:
		logging.info("Sub-directory {} already exists".format(sample_out_dir))

	#### create command line for emirge parameters #####
	logging.info("Preparing command for running EMIRGE with given parameters")
	parameter_cmd=""

	try:
		read1=open(read_1_file,"r")
		read1.close()
		read1_param="-1 {} ".format(read_1_file)
		parameter_cmd=parameter_cmd+read1_param
	except IOError:
		logging.error("Read1 file {} not found or cannot be opened".format(read_1_file))
		#print "Error: Read1 file not found or cannot be opened {}".format(read_1_file)
		sys.exit()

	try:
		read2=open(read_2_file,"r")
		read2.close()
		read2_param="-2 {} ".format(read_2_file)
		parameter_cmd=parameter_cmd+read2_param
	except IOError:
		logging.error("Read2 file {} not found or cannot be opened".format(read_2_file))
		#print "Error: Read2 file not found or cannot be opened {}".format(read_2_file)
		sys.exit()
		
	###### Prepare command parameter for EMIRGE #####
	logging.info("Checking parameters for EMIRGE")
	for parameter in emirge_parameter:

		if parameter=="PROGRAM":
			continue
		
		if parameter=="FASTA_DB":
			if emirge_parameter["FASTA_DB"]!=None:
				fasta_db_param="-f {} ".format(emirge_parameter["FASTA_DB"])
				parameter_cmd=parameter_cmd+fasta_db_param
				logging.info("--FASTA_DB/-f Ok")
			else:
				logging.error("Fasta database for EMIRGE not found at {}".format(emirge_parameter["FASTA_DB"]))
				#print "{} no fasta database found\n"
				sys.exit()

		if parameter=="BOWTIE_DB":
			if emirge_parameter["BOWTIE_DB"]!=None:
				bowtie_db_param="-b {} ".format(emirge_parameter["BOWTIE_DB"])
				parameter_cmd=parameter_cmd+bowtie_db_param
				logging.info("--BOWTIE_DB/-b Ok")
			else:
				logging.error("Bowtie database for EMIRGE not found at {}".format(emirge_parameter["BOWTIE_DB"]))
				#print "{} no bowtie database found\n"
				sys.exit()
				
		if parameter=="MAX_READ_LENGTH":
			if int(emirge_parameter["MAX_READ_LENGTH"]):
				max_read_len_param="-l {} ".format(emirge_parameter["MAX_READ_LENGTH"])
				parameter_cmd=parameter_cmd+max_read_len_param
				logging.info("--MAX_READ_LENGTH/-l Ok")
			else:
				logging.error("Invalid value {} used for max_read_length parameter".format(emirge_parameter["MAX_READ_LENGTH"]))
				#print "Invalid value for max_read_length parameter\n"
				sys.exit()
		

		if parameter=="INSERT_MEAN":
			if int(emirge_parameter["INSERT_MEAN"]):
				insert_mean_param="-i {} ".format(emirge_parameter["INSERT_MEAN"])
				parameter_cmd=parameter_cmd+insert_mean_param
				logging.info("--INSERT_MEAN/-i Ok")
			else:
				logging.error("Invalid value {} used for insert_mean parameter".format(emirge_parameter["INSERT_MEAN"]))
				#print "Invalid value for insert_mean parameter\n"
				sys.exit()
		

		if parameter=="INSERT_STDDEV":
			if int(emirge_parameter["INSERT_STDDEV"]):
				insert_stddev_param="-s {} ".format(emirge_parameter["INSERT_STDDEV"])
				parameter_cmd=parameter_cmd+insert_stddev_param
				logging.info("--INSERT_STDDEV/-s Ok")
			else:
				logging.error("Invalid value {} used for insert_stddev parameter".format(emirge_parameter["INSERT_STDDEV"]))
				#print "Invalid value for insert_stddev parameter\n"
				sys.exit()

		if parameter=="ITERATIONS":
			if int(emirge_parameter["ITERATIONS"]):
				iteration_param="-n {} ".format(emirge_parameter["ITERATIONS"])
				parameter_cmd=parameter_cmd+iteration_param
				logging.info("--ITERATIONS/-n Ok")
			else:
				logging.error("Invalid value {} used for number of iteration parameter".format(emirge_parameter["ITERATIONS"]))
				#print "Invalid value for number of iteration parameter\n"
				sys.exit()
			
		if parameter=="PROCESSORS":
			if int(emirge_parameter["PROCESSORS"]):
				if int(emirge_parameter["PROCESSORS"])!=cpu:
					logging.warning("Number of cpu cores {} specified for each job do not match with number of processors {} specified in the config file".format(cpu,emirge_parameter["PROCESSORS"]))
					#print "Number of cpu cores specified for each job do not match with number of processors specified in the config file\n"
					sys.exit()
				processor_param="-a {} ".format(emirge_parameter["PROCESSORS"])
				parameter_cmd=parameter_cmd+processor_param
				logging.info("--PROCESSORS/-a Ok")
			else:
				logging.error("Invalid value {} used for number of cpu cores".format(emirge_parameter["PROCESSORS"]))
				#print "Invalid value for number of cpu cores\n"
				sys.exit()

		if parameter=="MAPPING":
			try:
				bam_mapping_file=open(emirge_parameter["MAPPING"],"r")
				bam_mapping_file.close()
				bam_mapping_param="-m {} ".format(emirge_parameter["MAPPING"])
				parameter_cmd=parameter_cmd+bam_mapping_param
				logging.info("--MAPPING/-m Ok")
			except IOError:
				logging.error("BAM mapping file {} not exist or cannot be opened".format(emirge_parameter["MAPPING"]))
				#print "Error: BAM mapping file not exist or cannot be opened {}\n".format(emirge_parameter["MAPPING"])
				sys.exit()

		if parameter=="SNP_FRACTION_THRESH":
			if float(emirge_parameter["SNP_FRACTION_THRESH"]):
				snp_fraction_param="-p {} ".format(emirge_parameter["SNP_FRACTION_THRESH"])
				parameter_cmd=parameter_cmd+snp_fraction_param
				logging.info("--SNP_FRACTION_THRESH/-p Ok")
			else:
				logging.error("Invalid value {} used for SNP fraction threshold".format(emirge_parameter["SNP_FRACTION_THRESH"]))
				#print "Invalid value for SNP fraction threshold\n"
				sys.exit()

		if parameter=="VARIANT_FRACTION_THRESH":
			if float(emirge_parameter["VARIANT_FRACTION_THRESH"]):
				variant_fraction_param="-v {} ".format(emirge_parameter["VARIANT_FRACTION_THRESH"])
				parameter_cmd=parameter_cmd+variant_fraction_param
				logging.info("--VARIANT_FRACTION_THRESH/-v Ok")
			else:
				logging.error("Invalid value {} used for Variant fraction threshold".format(emirge_parameter["VARIANT_FRACTION_THRESH"]))
				#print "Invalid value for Variant fraction threshold\n"
				sys.exit()

		if parameter=="JOIN_THRESHOLD":
			if float(emirge_parameter["JOIN_THRESHOLD"]):
				rules=[float(emirge_parameter["JOIN_THRESHOLD"])>=0.95, float(emirge_parameter["JOIN_THRESHOLD"])<=1.00]
				if all(rules):
					join_threshold_param="-j {} ".format(emirge_parameter["JOIN_THRESHOLD"])
					parameter_cmd=parameter_cmd+join_threshold_param
					logging.info("--JOIN_THRESHOLD/-j Ok")
				else:
					logging.error("Value {} is out of range for join threshold".format(emirge_parameter["JOIN_THRESHOLD"]))
					#print "Value out of range for join threshold\n"
					sys.exit()
			else:
				logging.error("Invalid value {} used for join threshold parameter".format(emirge_parameter["JOIN_THRESHOLD"]))
				#print "Invalid value for join threshold parameter\n"
				sys.exit()

		if parameter=="MIN_DEPTH":
			if int(emirge_parameter["MIN_DEPTH"]):
				min_depth_param="-c {} ".format(emirge_parameter["MIN_DEPTH"])
				parameter_cmd=parameter_cmd+min_depth_param
				logging.info("--MIN_DEPTH/-c Ok")
			else:
				logging.error("Invalid value {} used for minimum depth".format(emirge_parameter["MIN_DEPTH"]))
				#print "Invalid value for minimum depth\n"
				sys.exit()

		if parameter=="NICE_MAPPING":
			if int(emirge_parameter["NICE_MAPPING"]):
				nice_mapping_param="--nice_mapping={} ".format(emirge_parameter["NICE_MAPPING"])
				parameter_cmd=parameter_cmd+nice_mapping_param
				logging.info("--NICE_MAPPING Ok")
			else:
				logging.error("Invalid input {} for nice mapping parameter".format(emirge_parameter["NICE_MAPPING"]))
				#print "Invalid input for nice mapping parameter\n"
				sys.exit()

		if parameter=="PHRED33":
			if emirge_parameter["PHRED33"].lower()=="yes":
				phred33_param="--phred33 "
				parameter_cmd=parameter_cmd+phred33_param
				logging.info("--PHRED33 Ok")
			else:
				logging.error("Invalid input {} for parameter --phred33".format(emirge_parameter["PHRED33"]))
				#print "Invalid input for parameter --phred33\n"
				sys.exit()

		if parameter=="SAVE_EVERY":
			if int(emirge_parameter["SAVE_EVERY"]):
				save_every_param="-e {} ".format(emirge_parameter["SAVE_EVERY"])
				parameter_cmd=parameter_cmd+save_every_param
				logging.info("--SAVE_EVERY/-e Ok")
			else:
				logging.error("Invalid input {} for --save_every parameter".format(emirge_parameter["SAVE_EVERY"]))
				#print "Invalid input for --save_every parameter\n"
				sys.exit()

		if parameter=="RESUME_FROM":
			if int(emirge_parameter["RESUME_FROM"]):
				resume_parameter="-r {} ".format(emirge_parameter["RESUME_FROM"])
				parameter_cmd=parameter_cmd+resume_parameter
				logging.info("--RESUME_FROM/-r Ok")
			else:
				logging.error("Invalid input {} for resume iteration parameter".format(emirge_parameter["RESUME_FROM"]))
				#print "Invalid input for resume iteration parameter\n"
				sys.exit()
	
	#### run Emirge command ####
	os.system("emirge.py {} {}".format(sample_out_dir, parameter_cmd))

	logging.info("Running EMIRGE for sample {}".format(sample_name))
	logging.info("emirge.py {} {}".format(sample_out_dir, parameter_cmd))

####### Function to rename emirge sequences for last iteration run #######
def rename_emirge_fasta(sample_name,project_path):

	logging.info("Analysis Step: Renaming/Analyzing EMIRGE-generated fasta sequence file")

	#### Get full path of the EMIRGE_RECONSTRUCT directory for give sample #####
	sample_emirge_dir=os.path.join(project_path,"EMIRGE_RECONSTRUCT")

	#### create output directory for the sample under processing ####
	sample_out_dir=os.path.join(project_path,"EMIRGE_RENAME")

	if not os.path.exists(sample_out_dir):
		os.makedirs(sample_out_dir)
		logging.info("Creating output sub-directory {}".format(sample_out_dir))
	else:
		logging.info("Sub-directory {} already exists".format(sample_out_dir))

	iteration_dir=get_dir_list(sample_emirge_dir)

	#### Iterate through list of emirge iteration directories in the sample output folder ###
	max_iter=00
	max_iter_dir=""

	logging.info("Finding directory for final iteration")
	
	for dir_name in iteration_dir:

		if dir_name=="initial_mapping":
			continue
		elif re.search("iter",dir_name):
			
			iteration=re.match(r"(iter)(\.)(\d+)",dir_name)
			
			iteration_number=iteration.group(3)

			if iteration_number>max_iter:

				max_iter=iteration_number
				max_iter_dir=dir_name
	

	logging.info("Number of iterations completed {}".format(max_iter))
	final_iter_dir=os.path.join(sample_emirge_dir,max_iter_dir)

	logging.info("Final iteration directory found at {}".format(final_iter_dir))

	renamed_fasta_file=os.path.join(sample_out_dir,sample_name+"_renamed.fasta")
	
	os.system("emirge_rename_fasta.py {} > {}".format(final_iter_dir, renamed_fasta_file))

	logging.info("Renaming fasta file in final iteration directory at {} to {}".format(final_iter_dir,renamed_fasta_file))
	logging.info("emirge_rename_fasta.py {} > {}".format(final_iter_dir, renamed_fasta_file))

	##### Reading EMIRGE reconstructed sequence file #####
	logging.info("Reading EMIRGE renamed sequence file at {}".format(renamed_fasta_file))
	records = list(SeqIO.parse(renamed_fasta_file, "fasta"))
	record_stats = [len(rec) for rec in records]

	emirge_reconstruction=os.path.join(sample_out_dir,"emirge_reconstruction.csv")
	logging.info("Logging information for EMIRGE-generated sequences in {}".format(emirge_reconstruction))
	emirge_seq_info=open(emirge_reconstruction,"w")
	emirge_seq_info.write("Sample_Name\tEmirge_Seq_Id\tSeq_Length\tPrior\tNormPrior\tSequence\n")
	#####Extracting information from EMIRGE reconstructed sequence file ####

	for rec in records:
		seq_id=rec.id
		seq_desc=rec.description
		seq=rec.seq
		seq_desc=seq_desc.rstrip('\n')
		seq_len=len(seq)
		

		parse_description=re.search(r"(\s)(Prior)(\=)(\d+\.\d+)(\s)",seq_desc)
		prior=parse_description.group(4)

		parse_description=re.search(r"(NormPrior)(\=)(.+)$",seq_desc)
		norm_prior=parse_description.group(3)

		emirge_seq_info.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(sample_name,seq_id,seq_len,prior,norm_prior,seq))

	emirge_seq_info.close()

	##### Summarizing results from EMIRGE Reconstruction ########
	summary_emirge_reconstruction=os.path.join(project_path,"summary_emirge_reconstruction.csv")
	logging.info("Preparing summary for EMIRGE-run in {}".format(summary_emirge_reconstruction))
	summary=open(summary_emirge_reconstruction,"w")

	summary.write("Sample_Name\tCount_Seq\tMax_Seq_Length\tMin_Seq_Length\tMean_Seq_Length\n")
	summary.write("{}\t{}\t{}\t{}\t{}\n".format(sample_name,len(records),max(record_stats),min(record_stats),statistics.mean(record_stats)))

	summary.close()


####### Function to perform nucleotide BLAST of predicted Emirge sequences against the reference database #####
def run_emirge_blast(sample_name,sample_dir,project_name,project_path,emirge_blast_parameters):

	logging.info("Analysis Step: Performing nucleotide BLAST for EMIRGE-generated sequences")
	blast_input_dir=os.path.join(sample_dir,project_name,"EMIRGE_RENAME")

	##### Print error message and terminate program if EMIRGE_RENAME directory is absent for given sample ####
	if not os.path.exists(blast_input_dir):
		logging.error("Directory {} does not exist or cannot be found".format(blast_input_dir))
		sys.exit()

	#### Read blast parameters from config dictonary #####
	blast_job_name=emirge_blast_parameters["BLAST_JOB_NAME"]
	blast_query=os.path.join(blast_input_dir,sample_name+"_renamed.fasta")
	blast_db=os.path.abspath(emirge_blast_parameters["DATABASE"])
	blast_program=emirge_blast_parameters["BLAST_PROGRAM"]
	blast_evalue=emirge_blast_parameters["EVALUE"]
	blast_identity=float(emirge_blast_parameters["PERCENT_IDENTITY_THRESHOLD"])
	blast_coverage=float(emirge_blast_parameters["QUERY_COVERAGE_THRESHOLD"])
	blast_target_num=emirge_blast_parameters["MAX_TARGET_SEQ"]
	blast_out_format=int(emirge_blast_parameters["OUTPUT_FORMAT"])
	blast_cpu=int(emirge_blast_parameters["NUM_THREADS"])
	only_parse_result=emirge_blast_parameters["ONLY_PARSE_RESULT"]
	add_taxonomy=emirge_blast_parameters["ADD_TAXONOMY"]
	taxon_dir=os.path.abspath(emirge_blast_parameters["TAXON_DIR"])

	logging.info("BLAST JOB NAME: {}".format(blast_job_name))
	logging.info("BLAST Query: {}".format(blast_query))
	logging.info("BLAST Database: {}".format(blast_db))

	##############   0   #  1 #  2 #  3   #  4 #  5 #  6   #  7   #  8 #  9   # 10 #  11  #  12  #  13    #  14  # 15  #  16  # 17     #  18    # 19 #  20  #  21  #  22   #   23    ######	
	blast_columns="qseqid qacc qlen sseqid sacc slen stitle qstart qend sstart send length evalue bitscore pident qcovs nident mismatch positive gaps qframe sframe staxids sskingdoms"

	db_name=os.path.basename(blast_db)

	#### create output directory for the sample under processing ####
	emirge_blast_out_dir=os.path.join(project_path,"EMIRGE_BLAST")
	blast_job_dir=os.path.join(emirge_blast_out_dir,blast_job_name)
	blast_raw_out_dir=os.path.join(emirge_blast_out_dir,blast_job_name,"RAW_BLAST")
	blast_parsed_out_dir=os.path.join(emirge_blast_out_dir,blast_job_name,"PARSED_BLAST")

	if not os.path.exists(emirge_blast_out_dir):
		os.makedirs(emirge_blast_out_dir)
		logging.info("Creating output sub-directory EMIRGE_BLAST in project directory {}".format(project_path))
	else:
		logging.info("Output sub-directory EMIRGE_BLAST in project directory {} already exist".format(project_path))

	if not os.path.exists(blast_job_dir):
		os.makedirs(blast_job_dir)
		logging.info("Creating BLAST job directory {}".format(blast_job_dir))
	else:
		logging.info("BLAST job directory {} already exist".format(blast_job_dir))

	if not os.path.exists(blast_raw_out_dir):
		os.makedirs(blast_raw_out_dir)
		logging.info("Creating BLAST output directory for storing raw result {}".format(blast_raw_out_dir))
	else:
		logging.info("BLAST output directory for storing raw result {} already exist".format(blast_raw_out_dir))

	if not os.path.exists(blast_parsed_out_dir):
		os.makedirs(blast_parsed_out_dir)
		logging.info("Creating BLAST output directory for storing parsed result {}".format(blast_parsed_out_dir))
	else:
		logging.info("BLAST output directory for storing parsed result {} already exist".format(blast_parsed_out_dir))

	blast_raw_output_file=os.path.join(blast_raw_out_dir,sample_name+"_blast_table.tsv")
	blast_parsed_output_file=os.path.join(blast_parsed_out_dir,sample_name+"_blast_parsed_table.tsv")

	#### Execute BLAST Command on the sample query ######
	if only_parse_result.upper()=="NO":

		logging.info("Running nucleotide BLAST Search for {} against {} database".format(blast_query,blast_db))

		os.system("{} -query {} -db {} -evalue {} -max_target_seqs {} -num_threads {} -out {} -outfmt '{} {}'".format(blast_program,
                                           					                                  	blast_query,
					                                                                          	blast_db,
                                                                                                                  	blast_evalue,
                                                                                                                 	blast_target_num,
                                                                                                                  	blast_cpu,
                                                                                                                  	blast_raw_output_file,
                                                                                                                 	blast_out_format,
                                                                                                                  	blast_columns))

		logging.info("{} -query {} -db {} -evalue {} -max_target_seqs {} -num_threads {} -out {} -outfmt '{} {}'".format(blast_program,
                                           					                                  	blast_query,
					                                                                          	blast_db,
                                                                                                                  	blast_evalue,
                                                                                                                  	blast_target_num,
                                                                                                                 	blast_cpu,
                                                                                                                 	blast_raw_output_file,
                                                                                                                 	blast_out_format,
                                                                                                                 	blast_columns))

		logging.info("Location of raw BLAST output file for sample {} is {}".format(sample_name,blast_raw_output_file))
	
	##### Parse BLAST RESULTS ######
	logging.info("Parsing BLAST result from {}".format(blast_raw_output_file))

	os.system("python runEmirgeParseBlast.py --sample_name {} --project_name {} --project_path {} --blast_file {} --db_name {} --emirge_file {} --output_file {} --identity_cutoff {} --qcov_cutoff {} --add_taxonomy {} --taxon_dir {}".format(sample_name, 
																														    project_name, 
																														    project_path,
																														    blast_raw_output_file,
																			      											    blast_db,
																			      											    blast_query,
																			      											    blast_parsed_output_file,
																			      											    blast_identity,
																			      											    blast_coverage,
																			      											    add_taxonomy,
																			      											    taxon_dir))

	logging.info("python runEmirgeParseBlast.py --sample_name {} --project_name {} --project_path {} --blast_file {} --db_name {} --emirge_file {} --output_file {} --identity_cutoff {} --qcov_cutoff {} --add_taxonomy {} --taxon_dir {}".format(sample_name, 
																														    project_name, 
																														    project_path,
																														    blast_raw_output_file,
																			      											    blast_db,
																			      											    blast_query,
																			      											    blast_parsed_output_file,
																			      											    blast_identity,
																			      											    blast_coverage,
																			      											    add_taxonomy,
																			      											    taxon_dir))
	logging.info("Parsed BLAST output file for sample {} is at {}".format(sample_name,blast_parsed_output_file))

##### Function to select best blast hit from priority database list #####

def run_emirge_select(sample_name,sample_dir,project_name,project_path,emirge_select_parameters):


	logging.info("Analysis Step: Selecting best-hit from the parsed BLAST result in order of database priority")

	input_dir=os.path.join(project_path,"EMIRGE_BLAST")
	emirge_select_out_dir=os.path.join(project_path,"EMIRGE_SELECT")
	emirge_select_out_file=os.path.join(emirge_select_out_dir,sample_name+"_blast_besthit_table.tsv")
	priority_db_list=emirge_select_parameters["PRIORITY_DB_REF"]
	num_hit=emirge_select_parameters["NUM_HIT_PER_QUERY"]
	
	logging.info("Order of database priority: {}".format(priority_db_list))

	if not os.path.exists(emirge_select_out_dir):
		os.makedirs(emirge_select_out_dir)
		logging.info("Create output sub-directory for selected blast-hit at {}".format(emirge_select_out_dir))
	else:
		logging.info("Output sub-directory for selected blast-hit at {} already exist".format(emirge_select_out_dir))

	
	logging.info("Running selection of best BLAST hit from {}".format(input_dir))

	os.system("python runEmirgeSelectBlastHit.py --sample_name {} --project_name {} --project_path {} --input_dir {} --output_file {} --priority_db_ref {} --num_hits {}".format(sample_name,
																						     project_name,
																						     project_path,
																						     input_dir,
															       							     emirge_select_out_file,
															       							     priority_db_list,
															       							     num_hit))

	logging.info("python runEmirgeSelectBlastHit.py --sample_name {} --project_name {} --project_path {} --input_dir {} --output_file {} --priority_db_ref {} --num_hits {}".format(sample_name,
																						     project_name,
																						     project_path,
																						     input_dir,
															       							     emirge_select_out_file,
															       							     priority_db_list,
															       							     num_hit))

	logging.info("Best BLAST hit output file for sample {} is at {}".format(sample_name,emirge_select_out_file))

			
			
if __name__ == '__main__':
	main()
	
	






