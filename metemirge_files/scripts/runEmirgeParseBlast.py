#!/usr/bin/env python

"""
Shalabh Thakur
May, 2018
Sycuro Lab, University of Calgary
"""
"""
Description: This script is used to parse BLAST results for EMIRGE sequences.
"""
import sys
import os
import glob
import argparse
import csv
import re
import shutil
import logging
import datetime
import subprocess
import copy
import sqlite3
from argparse import RawTextHelpFormatter
from collections import defaultdict
from Bio import SeqIO
from ete3 import NCBITaxa



##### MAIN FUNCTION #######
def main():

	"""
	Main Function
	"""
	now = datetime.datetime.now()
	
	c_args = parse_args(__file__)

	sample_name=c_args['sample_name']
	project_name=c_args['project_name']
	project_path=os.path.abspath(c_args['project_path'])
	blast_file=os.path.abspath(c_args['blast_file'])
	db_name=c_args['db_name']
	emirge_file=os.path.abspath(c_args['emirge_file'])
	output_file=os.path.abspath(c_args['output_file'])
	identity_cutoff=c_args['identity_cutoff']
	qcov_cutoff=c_args['qcov_cutoff']
	add_taxonomy=c_args['add_taxonomy']
	taxonomy_dir=(os.path.abspath(c_args['taxon_dir']))
	silva_taxonomy_db=c_args['silva_taxadb']
	ncbi_taxonomy_db=c_args['ncbi_taxadb']

	log_file=os.path.join(project_path,sample_name+"_"+project_name.lower()+"_run.log")
	logging.basicConfig(filename=log_file,format="[%(asctime)s] %(levelname)-8s %(message)s",datefmt='%a, %d %b %Y %H:%M:%S',level=logging.INFO)

	#### Download NCBI accession2taxid mapping files for the first time script is run,
	#### If the accession2taxid mapping files from ncbi already present then skip downloading
	#os.system("python update_ncbi_accn2taxid.py --db_dir {} --create".format(os.path.abspath(c_args['taxon_dir'])))

	####Calculate normalized prior and emirge sequence #####
	emirge_seq_record=get_norm_prior_seq(emirge_file)

	#### Calculate full length percent identity and parse blast resut ####
	formatted_blast_table=parse_blast(blast_file,emirge_seq_record,db_name)

	####Fix multiple HSPs #####
	fixed_hsp_blast_table=fix_multi_hsp(formatted_blast_table)

	####Convert HSPs to HIT ######
	formatted_hit_blast_table=convert_hsp2hit(fixed_hsp_blast_table,identity_cutoff,qcov_cutoff,add_taxonomy,taxonomy_dir,silva_taxonomy_db,ncbi_taxonomy_db)

	#### Write Formatted Blast Table to the output file ####
	write_blast_table(formatted_hit_blast_table,output_file)

	

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
	parser.add_argument("--sample_name",help="name of the sample", required=True)
	parser.add_argument("--project_name",help="name of the project", required=True)
	parser.add_argument("--project_path",help="path of the project", required=True)
	parser.add_argument("--blast_file", help="name and path of the blast result file", required=True)
	parser.add_argument("--db_name", help="name of blast database used for EMIRGE", required=True)
	parser.add_argument("--emirge_file", help="name and path of renamed emrige fasta sequence file", required=True)
	parser.add_argument("--output_file", help="name path to the parsed output file", required=True)
	parser.add_argument("--identity_cutoff", type=float, default=20.0, help="cut-off for percent identity in the blast result")
	parser.add_argument("--qcov_cutoff", type=float, default=20.0, help="cut-off for percent query coverage in the blast result")
	parser.add_argument("--add_taxonomy",default="yes", help="specify if taxonomic information should be added to blast result")
	parser.add_argument("--taxon_dir", default="../taxonomy",help="name and path of NCBI taxonomy dir")
	parser.add_argument("--silva_taxadb", default="silva_taxonomy_v132.sqlite", help="name of the SILVA taxonomy sqlite database")
	parser.add_argument("--ncbi_taxadb", default="ncbi_taxonomy.sqlite", help="name of the NCBI taxonomy sqlite database")
	
	args=parser.parse_args()

	return vars(args)



def parse_blast(blast_file, emirge_seq_record, db_name):

	"""
	Reads in tabular BLAST results, calculate percent identity for full length query sequence and create a dictionary for each row and column

	param blast_file: full path and name of blast output file

	return blast_table: dict of blast columns sorted by query_id/subject_id pair
		
	The method/formula of calculating full length pident was taken from https://github.com/DenefLab/EMIRGE/blob/master/calc_full_length_pident.R
	"""

	#### Read Blast Raw output file #####
	blast_result=open(blast_file,"r")

	formatted_blast_table=defaultdict(list)

	#print "[PARSER_STEP 2]: Reading BLAST results from {}\n".format(blast_file)

	logging.info("[PARSER_STEP 2]: Reading BLAST results from {}".format(blast_file))

	for index, blast_row in enumerate(blast_result):

		blast_row=blast_row.rstrip("\n")
		blast_column=blast_row.split("\t")

		query_subject_pair="{}:{}".format(blast_column[0],blast_column[3])
		q_norm_prior=float(emirge_seq_record[blast_column[0]]['norm_prior'])
		q_seq=emirge_seq_record[blast_column[0]]['sequence']

		#### calculate full length percent identity ####
		#full_len_pid=calc_full_pid(blast_column)

		#blast_column.append(full_len_pid)
		blast_column.append(q_norm_prior)
		blast_column.append(db_name)
		blast_column.append(q_seq)

		hsp_hit={}

		hsp_hit['qseqid'] = blast_column[0]
		hsp_hit['qacc'] = blast_column[1]
	 	hsp_hit['qlen'] = int(blast_column[2])
		hsp_hit['sseqid'] = blast_column[3]
		hsp_hit['sacc'] = blast_column[4]
		hsp_hit['slen'] = int(blast_column[5])
		hsp_hit['stitle'] = blast_column[6]
		hsp_hit['qstart'] = int(blast_column[7])
		hsp_hit['qend'] = int(blast_column[8])
		hsp_hit['sstart'] = int(blast_column[9])
		hsp_hit['send'] = int(blast_column[10])
		hsp_hit['length'] = int(blast_column[11])
		hsp_hit['evalue'] = blast_column[12]
		hsp_hit['bitscore'] = float(blast_column[13])
		hsp_hit['pident'] = float(blast_column[14])
		hsp_hit['qcovs'] = float(blast_column[15])
		hsp_hit['nident'] = int(blast_column[16])
		hsp_hit['mismatch'] = int(blast_column[17])
		hsp_hit['positive'] = int(blast_column[18])
		hsp_hit['gaps'] = int(blast_column[19])
		hsp_hit['qframe'] = blast_column[20]
		hsp_hit['sframe'] = blast_column[21]
		hsp_hit['staxid'] = blast_column[22]
		hsp_hit['skingdom'] = blast_column[23]
		hsp_hit['q_norm_prior'] = blast_column[24]
		hsp_hit['db_name'] = blast_column[25]
		hsp_hit['q_seq'] = blast_column[26]
		hsp_hit['is_species'] = False
		hsp_hit['taxonomy'] = None

		formatted_blast_table[query_subject_pair].append(hsp_hit)
			
	return(formatted_blast_table)


def calc_full_pid(blast_column):

	pident=float(blast_column[14])
	aln_length=int(blast_column[11])
	q_align=(int(blast_column[8]) - int(blast_column[7]))+1
	qlen=int(blast_column[2])

	full_len_pid=round((pident * aln_length)/ (aln_length - q_align + qlen),2)

	return(full_len_pid)


def get_norm_prior_seq(emirge_file):

	"""
	Reads header for each sequence from renamed emirge fasta file and extract normalized prior/relative abundance
	Also add sequence for corresponding subject from the reference database

	param emirge_file: name and path of renamed emirge fasta file
	param formatted_blast_table: dictonary of blast result table with full length pid

	return formatted_blast_table added with normalized prior and subject sequence
	"""

	emirge_seq_record={}
	#### read emirge fasta file ####

	#print "[PARSER_STEP 1]: Getting Normalize Prior and Sequence information from EMIRGE file\n"
	
	logging.info("[PARSER_STEP 1]: Getting Normalize Prior and Sequence information from {}".format(emirge_file))

	with open(emirge_file, "r") as fasta_handle:
        	seq_record = SeqIO.parse(fasta_handle, "fasta")

		for record in seq_record:
			
			q_id=record.id
			q_seq=record.seq
			q_desc=record.description

			q_desc=q_desc.rstrip('\n')

			parse_description=re.search(r"(NormPrior)(\=)(.+)$",q_desc)
			q_norm_prior=parse_description.group(3)
			
			emirge_seq_record[q_id]={'norm_prior':q_norm_prior, 'sequence':q_seq}

	return(emirge_seq_record)


def fix_multi_hsp(formatted_blast_table):

	"""
	Reads formatted blast result, find compatible hsps and re-calculate identities and bitscore for multiple hsp hits
	"""
	fixed_hsp_blast_table={}

	logging.info("[PARSER_STEP 3]: Fixing multi-hsp hits")

	for query_subject_pair in formatted_blast_table:

		#sorted_hsp_list=sorted(formatted_blast_table[query_subject_pair], key=lambda k: k['qstart'], reverse=False)
	
		##### find compatible hsps with highest bitscore that do not cross and overlap each other ##
		compatible_hsp=list()

		logging.info("[PARSER_STEP 3]: Computing best hsp for {}".format(query_subject_pair))

		compatible_hsp = compute_best_hsp(formatted_blast_table[query_subject_pair],hsp_list=list())	

		fixed_hsp_blast_table[query_subject_pair]=compatible_hsp

	return(fixed_hsp_blast_table)


# Return summed bitscores of best combo. To figure out which HSPs
# correspond to this optimal score, we must also store whether the optimal
# solution includes the last HSP in `hsps` at each step, then recursively
# reconstruct the optimal solution after the algorithm finishes. See the
# linked code for details.

def compute_best_hsp(hsps,hsp_list=list()):
  	"""
  	Memoize the function so that, if we have already computed the solution
  	for `hsps`, just fetch and return that value immediately without
  	performing below computations.

	Based on code found here:
    	http://jeff.wintersinger.org/posts/2014/07/designing-an-algorithm-to-
    	compute-the-optimal-set-of-blast-hits/
  	"""
	if len(hsps) == 0:
		return hsp_list
	if len(hsps) == 1:
    		# Trivial solution: one HSP, so optimal solution is just itself.
		hsp_list.append(hsps[0])
		return hsp_list

  	# Last HSP
  	last_hsp = hsps[-1]
  	# All HSPs except last
  	previous = hsps[:-1]

  	# Find subset of HSPs in `previous` that don't overlap `last_hsp`.
  	compatible = find_compatible(last_hsp, previous)

	without_list = copy.deepcopy(hsp_list)
    	with_list = copy.deepcopy(hsp_list)
    	with_list.append(last_hsp)

  	best_without_last = compute_best_hsp(previous,without_list)
	best_with_last    = compute_best_hsp(compatible, with_list)

	logging.info("[PARSER_STEP 3]: Finding compatible hsps from multi-hsp hits")
	compatible_hsp=max([best_without_last, best_with_last],key=lambda (x): sum([float(y['bitscore']) for y in x]))
	
	return compatible_hsp

##### Find Compatible HSP #####
def find_compatible(target, hsps):
  	'''Find hsps amongst `hsps` that don't overlap or cross `target`. They
  	may not be mutually compatible, as they are only guaranteed to be
  	compatible with `target`.'''
	compatible = []

  	for hsp in hsps:
    	# Don't define target as being compatible with itself.
    		if hsp == target:
      			continue

   		if target['qstart'] <= hsp['qstart']:
      			first, second = target, hsp
    		else:
      			first, second = hsp, target

    		overlap = (second['qstart'] <= first['qend'] or
              		second['sstart'] <= first['send'])
    		if not overlap:
      			compatible.append(hsp)

	return compatible

##### convert HSP to HIT format ####

def convert_hsp2hit(fixed_hsp_blast_table,identity_cutoff,qcov_cutoff,add_taxonomy,taxonomy_dir,silva_taxonomy_db,ncbi_taxonomy_db):

	formatted_hit_blast_table={}

	logging.info("[PARSER_STEP 4]: Converting HSP format into Hit format")

	for query_subject_pair in fixed_hsp_blast_table:

		multi_hsp_hit=list()
		hsp_hit_list=fixed_hsp_blast_table[query_subject_pair]
	
		total_subject_pident=hsp_hit_list[0]['pident']
		total_subject_qcov=hsp_hit_list[0]['qcovs']
		total_subject_bitscore=hsp_hit_list[0]['bitscore']
		total_query_len=hsp_hit_list[0]['qlen']
		is_species=False
		is_multi_hsp=False
		num_hsp=len(hsp_hit_list)
		subject_accession=hsp_hit_list[0]['sacc']
		taxid=hsp_hit_list[0]['staxid']
		stitle=hsp_hit_list[0]['stitle']
		taxonomy_path="unclassified"
		taxonomy_source="unknown"
		ncbi_taxid="N/A"
		silva_taxid="N/A"
		rank="N/A"
		subject_species="unclassified"
		genus="unclassified"


		####Clean Accession id #####
		if re.search(r"^(ref)(\|)",subject_accession):
			accession_line=re.search(r"^(\w+)(\|)(\w+)(\|)(\w+)(\|)(\w+)(\|)",subject_accession)
			subject_accession=accession_line.group(3)

		#### Get Name of Subject Species from Cleaned Sequence Title #######
		stitle=stitle.replace(hsp_hit_list[0]['sseqid'],"")
		stitle=re.sub(r"^\s+","",stitle)

		species_description=re.search(r"^(\S+)(\s)(\S+)",stitle)

		if re.search(r"[\,\.\?\+\-\!\:\_\;]",species_description.group(1)):
			subject_species=species_description.group(1)
			subject_species=re.sub(r"[\,\.\?\+\-\!\:\;]","",subject_species)
			subject_species=re.sub(r"\_"," ",subject_species)
		else:
			subject_species=species_description.group(1)+" "+species_description.group(3)

		##### Get Name of the Genus ####
		subject_species=re.sub(r"[\[\]\(\)]","",subject_species)

		print "SPECIES:",subject_species

		genus_description=re.search(r"^(\w+)",subject_species)
		genus=genus_description.group(1)

		print "GENUS:",genus

		##### search for taxonomic information for the hit ######
		##### if taxid is assigned during blast search 
		##### then directory search for taxonomy classfication using ETE2 taxonomy
		##### based on taxid. However, if no taxid assigned, then use genus/species/strain 
		##### names from stitle to search for taxonomy

		if add_taxonomy=="yes":
			taxonomy_path,taxonomy_source,ncbi_taxid,silva_taxid,rank=get_taxonomy(genus,subject_accession,taxid,stitle,taxonomy_dir,silva_taxonomy_db,ncbi_taxonomy_db)

		print taxonomy_path

		if len(hsp_hit_list)>1:

			sum_num_ident=sum(k['nident'] for k in hsp_hit_list)
			sum_aln_length=sum(k['length'] for k in hsp_hit_list)
			total_subject_bitscore=sum(k['bitscore'] for k in hsp_hit_list)
			total_subject_pident=round((int(sum_num_ident) / int(sum_aln_length))*100,2)
			total_subject_qcov=round((int(sum_aln_length) / int(total_query_len))*100,2)
			is_multi_hsp=True


		if float(total_subject_pident)>=float(identity_cutoff) and float(total_subject_qcov)>=float(qcov_cutoff):
			is_species=True

		multi_hsp_hit.append(hsp_hit_list[0]['qseqid'])
		multi_hsp_hit.append(hsp_hit_list[0]['qacc'])
		multi_hsp_hit.append(hsp_hit_list[0]['qlen'])
		multi_hsp_hit.append(hsp_hit_list[0]['sseqid'])
		multi_hsp_hit.append(subject_accession)
		multi_hsp_hit.append(hsp_hit_list[0]['slen'])
		multi_hsp_hit.append(hsp_hit_list[0]['stitle'])
		multi_hsp_hit.append(subject_species)
		multi_hsp_hit.append(hsp_hit_list[0]['q_norm_prior'])
		multi_hsp_hit.append(hsp_hit_list[0]['db_name'])
		multi_hsp_hit.append(taxonomy_path)
		multi_hsp_hit.append(ncbi_taxid)
		multi_hsp_hit.append(silva_taxid)
		multi_hsp_hit.append(rank)
		multi_hsp_hit.append(taxonomy_source)
		multi_hsp_hit.append(is_species)
		multi_hsp_hit.append(is_multi_hsp)
		multi_hsp_hit.append(num_hsp)
		multi_hsp_hit.append('_'.join(map(str,(hsp_hit['qstart'] for hsp_hit in  hsp_hit_list))))
		multi_hsp_hit.append('_'.join(map(str,(hsp_hit['qend'] for hsp_hit in  hsp_hit_list))))
		multi_hsp_hit.append('_'.join(map(str,(hsp_hit['sstart'] for hsp_hit in  hsp_hit_list))))
		multi_hsp_hit.append('_'.join(map(str,(hsp_hit['send'] for hsp_hit in  hsp_hit_list))))
		multi_hsp_hit.append('_'.join(map(str,(hsp_hit['length'] for hsp_hit in  hsp_hit_list))))
		multi_hsp_hit.append('_'.join(map(str,(hsp_hit['evalue'] for hsp_hit in  hsp_hit_list))))
		multi_hsp_hit.append(total_subject_bitscore)
		multi_hsp_hit.append('_'.join(map(str,(hsp_hit['bitscore'] for hsp_hit in  hsp_hit_list))))
		multi_hsp_hit.append(total_subject_pident)
		multi_hsp_hit.append('_'.join(map(str,(hsp_hit['pident'] for hsp_hit in  hsp_hit_list))))
		multi_hsp_hit.append(total_subject_qcov)
		multi_hsp_hit.append('_'.join(map(str,(hsp_hit['qcovs'] for hsp_hit in  hsp_hit_list))))
		multi_hsp_hit.append('_'.join(map(str,(hsp_hit['nident'] for hsp_hit in  hsp_hit_list))))
		multi_hsp_hit.append('_'.join(map(str,(hsp_hit['mismatch'] for hsp_hit in  hsp_hit_list))))
		multi_hsp_hit.append('_'.join(map(str,(hsp_hit['positive'] for hsp_hit in  hsp_hit_list))))
		multi_hsp_hit.append('_'.join(map(str,(hsp_hit['gaps'] for hsp_hit in  hsp_hit_list))))
		multi_hsp_hit.append('_'.join(map(str,(hsp_hit['qframe'] for hsp_hit in  hsp_hit_list))))
		multi_hsp_hit.append('_'.join(map(str,(hsp_hit['sframe'] for hsp_hit in  hsp_hit_list))))
		multi_hsp_hit.append(hsp_hit_list[0]['q_seq'])

		formatted_hit_blast_table[query_subject_pair]=multi_hsp_hit	

	return(formatted_hit_blast_table)


##### get taxonomy function ######

def get_taxonomy(genus,accession,taxid,stitle,taxonomy_dir,silva_taxonomy_db,ncbi_taxonomy_db):

	"""
	get the taxonomy information based on taxid or genus name
	"""

	logging.info("[PARSER_STEP 5]: Assigning Taxonomy for {}".format(accession))

	ncbi= NCBITaxa()

	taxonomy_path="unclassified"
	taxonomy_source="unknown"
	rank="N/A"
	silva_taxid="N/A"
	ncbi_taxid="N/A"

	silva_taxonomy=os.path.join(taxonomy_dir,silva_taxonomy_db)
	ncbi_taxonomy=os.path.join(taxonomy_dir,ncbi_taxonomy_db)
	#ncbi_accn2taxid=os.path.join(taxonomy_dir,"ncbi_accn2taxid")

	conn_silva_db=sqlite3.connect(silva_taxonomy)
	conn_ncbi_db=sqlite3.connect(ncbi_taxonomy)
	
	### if taxid is not assigned during blast, then use accession number to find
	### taxid from ncbi accn2taxid mapping files 
	accession_number=accession

	if re.search(r"^(ref)(\|)",accession):
		accession_line=re.search(r"^(\w+)(\|)(\w+)(\|)(\w+)(\|)(\w+)(\|)",accession)
		accession_number=accession_line.group(3)

	#### remove version information from accession number ####
	accession_number=re.sub("\..+","",accession_number)

	db_record_tuple=(accession_number,)

	###### First check taxonomy in Silva taxonomy ######
	###### If sequences is annotated in Silva give it a first priority ####
	silva_db_cursor=conn_silva_db.cursor()

	silva_db_cursor.execute("SELECT distinct(taxa_path),taxid,db_version from taxonomy where accession_id=?",db_record_tuple)
	taxonomy_line=silva_db_cursor.fetchone()

	if taxonomy_line:
		#### If taxonomy found in SILVA use it ####
		taxonomy_path=taxonomy_line[0]
		silva_taxid=taxonomy_line[1]
		silva_version=taxonomy_line[2]
		taxonomy_source="Silva_db_v{}".format(silva_version)
	else:
		#### Check in NCBI taxonomy when there is no match in SILVA database ####
		#cmd_ncbi=["grep","-m1","-w","-r",'{}'.format(accession_number),ncbi_accn2taxid]
		#taxonomy_ncbi=subprocess.Popen(cmd_ncbi,stdout=subprocess.PIPE)
		#taxonomy_line=taxonomy_ncbi.communicate()[0]

		ncbi_db_cursor=conn_ncbi_db.cursor()

		ncbi_db_cursor.execute("SELECT distinct(ncbi_taxid),ncbi_accession_id from accession2taxid where ncbi_accession_id=?",db_record_tuple)
		taxid_line=ncbi_db_cursor.fetchone()

		if taxid_line:
			taxid = taxid_line[0]
			taxonomy_source="Ncbi_db"
	
		if not (taxid=="N/A" or taxid==None):

			taxonomy_path=get_ncbi_taxonomy(taxid)

		else:

			if not re.search("Unculture",genus,re.IGNORECASE):

				genus_taxid=ncbi.get_name_translator([genus])

				if genus_taxid:
					taxid=genus_taxid[genus].pop()
					taxonomy_path=get_ncbi_taxonomy(taxid)
					taxonomy_source="Ncbi_db"

			else:	
				taxonomy_path="unclassified"

	taxonomy_path=re.sub(r"^;","",taxonomy_path)
	taxonomy_path=re.sub(r";$","",taxonomy_path)

	##### get rank for final taxonomic level ######
	taxonomic_column=taxonomy_path.split(";")
	last_tax_level=taxonomic_column.pop()
	rank,ncbi_taxid=get_ncbi_taxa_rank(last_tax_level)

	print taxonomy_path
	
	return(taxonomy_path,taxonomy_source,ncbi_taxid,silva_taxid,rank)


##### get ncbi taxonomy for given taxonomic Id ######

def get_ncbi_taxonomy(taxid):

	ncbi= NCBITaxa()

	lineage = ncbi.get_lineage(taxid)
	names = ncbi.get_taxid_translator(lineage)
	ranks = ncbi.get_rank(lineage)
	ncbi_taxonomy_path=""

	for taxid in lineage:

		if not ranks[taxid]=="no rank":
			ncbi_taxonomy_path = ncbi_taxonomy_path +";"+names[taxid]

	return(ncbi_taxonomy_path)

#### get ncbi rank for taxonomic name ####

def get_ncbi_taxa_rank(taxa_name):

	ncbi= NCBITaxa()

	name2taxid=ncbi.get_name_translator([taxa_name])
	rank="N/A"
	ncbi_taxid="N/A"

	if name2taxid:
		ncbi_taxid=name2taxid[taxa_name].pop()
		ncbi_ranks=ncbi.get_rank([ncbi_taxid])
		rank=ncbi_ranks[ncbi_taxid]
	
	return(rank,ncbi_taxid)

#### Write Parsed Blast Table #####

def write_blast_table(formatted_hit_blast_table,output_file):

	"""
	Reads formatted blast result dictonary and write it to output file
	"""	
	logging.info("[PARSER_STEP 6]: Writing parsed BLAST result in {}".format(output_file))

	blast_report=open(output_file,"w")

	##### Write Column header ######

	#print "[PARSER_STEP 3]: Writing Formatted Blast Results to the output file {}\n".format(output_file)

	blast_report.write("QUERY_ID\tQUERY_ACCESSION\tQUERY_LENGTH\tSUBJECT_ID\tSUBJECT_ACCESSION\tSUBJECT_LENGTH\tSUBJECT_TITLE\tSUBJECT_SPECIES\tEMIRGE_NORM_PRIOR\tBLAST_DB\tTAXONOMY\tNCBI_TAXID\tSILVA_TAXID\tRANK\tTAXONOMY_SOURCE\tIS_SPECIES\tIS_MULTI_HSP\tNUM_HSP\tQUERY_START\tQUERY_END\tSUBJECT_START\tSUBJECT_END\t"
			   "ALIGNMENT_LENGTH\tEVALUE\tTOTAL_BITSCORE\tHSP_BITSCORE\tTOTAL_PERCENT_IDENTITY\tHSP_PERCENT_IDENTITY\tTOTAL_QUERY_COVERAGE\tHSP_QUERY_COVERAGE\tNUMBER_IDENTITIES\tMISMATCHES\tPOSITIVE\tGAPS\tQUERY_FRAME\t"
			   "SUBJECT_FRAME\tSEQUENCE\n")

	for query_subject_pair in formatted_hit_blast_table:

			hit_column=formatted_hit_blast_table[query_subject_pair]

			column_value="\t".join(map(str,(hit_column)))

			#print column_value

			blast_report.write("{}\n".format(column_value))

		

if __name__ == '__main__':
	main()
