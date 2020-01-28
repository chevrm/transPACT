## Copyright (c) 2016 Xiaowen Lu
## Wageningen University
## Bioinformatics Group

import re




from antismash.utils import sort_XonY

#--------------------------------
#-- function: get_strand_info
#--------------------------------
def get_strand_info(domain_info):
	"""
	This the the function to get the strand info of KS domains on a gene clusters. the strand info are order according to
	the DNA location of the gene. And the second column contains information for whether a TE domain is identified on the BGC
	:param file_dir: a string, which is the directory name of the folder containing KS domain of all BGCs
	:param output: a string, which is the file name of the output
	:return: NA. the output is written to a txt file

	file_dir = "/Users/Xiaowen/Documents/trans_AT/KS_domains_PredtransAT_hmm2_20160330"
	output = "/Users/Xiaowen/Documents/trans_AT/KS_domains_PredtransAT_output_20160330/KS_domains_strand_info.txt"
	get_strand_info(file_dir, output)
	"""
	gene_loc = [int(a.split("|")[2].split("-")[0]) for a in domain_info]
	strand = [a.split("|")[3] for a in domain_info]
	strand_dict = {}
	for x, y in zip(strand, gene_loc):
		strand_dict.setdefault(y, []).append(x)

	key_sorted = sorted(strand_dict.keys())
	strand_sorted = []
	for k in key_sorted:
		strand_sorted = strand_sorted + strand_dict[k]

	strand_info = ''.join(strand_sorted)

	return(strand_info)



#--------------------------------
#-- IndexDomain
#--------------------------------
def index_domain(accession, strand_info, domain_ab):
	'''
	This is the function to order the KS domains in an assembly line. Since the KS domain can be encoded either on a positive or negative strand, which has an effect
	for the ordering of the KS domains. Here, we order the KS according to the following rules:
	i) if all KS domains are encoded in only positive strand or only negative strand. then we order the KS domain in the same order of the strand location where the KS is in.
	(Please refer to note for details https://www.evernote.com/shard/s311/nl/2147483647/05fabd75-6ac2-48c4-b916-95e388ee0bc9/)
	ii) if the KS domains are partially encoded on negative strand. a) If the KS domains are encoded as ----+++ or +++---, the KS domains on '+' are ordered as 1,2,... and the ones
	on '-' are order as 1b,2b,... At this step, we also make use of the location of TE domain to determine the KS domain order.
	b) If the KS domains are encoded as mixture ----+---- or +++-+++, no index is assigned.
	:param input:
	:param output:
	:param domain_ab:
	:param KS_strand_info:
	:return:

	input = "/Users/Xiaowen/Documents/trans_AT/KS_domains_PredtransAT_hmm2_20160330/AMWE01000004_c3_KS_indexed.txt"
	output = "/Users/Xiaowen/Documents/trans_AT/KS_domains_PredtransAT_hmm2_index_20160330/AMWE01000004_c3_KS_indexed.txt"
	KS_strand_info = "/Users/Xiaowen/Documents/trans_AT/KS_domains_PredtransAT_output_20160330/KS_domains_strand_info.txt"
	TE_folder = "/Users/Xiaowen/Documents/trans_AT/TE_domains_PredtransAT_hmm2"
	domain_ab = "KS"
	'''

	accession_common = ["|".join(a.split("|")[0:3]) for a in accession]
	accessiondict = {}
	for x,y in zip(accession, accession_common):
		accessiondict.setdefault(y.split('|')[2].split('-')[0],[]).append(x)

	strand_info_bgc_neg = strand_info.count('-')
	strand_info_bgc_pos = strand_info.count('+')

	if strand_info_bgc_neg == 0:
		domain_name_dict = {}
		keys_sorted = sorted([int(k) for k in accessiondict.keys()])
		domnr = 1
		for key in keys_sorted:
			domain_index = [k.split('|')[-1].replace(domain_ab,'') for k in accessiondict.get(str(key))]
			val = sort_XonY(accessiondict.get(str(key)), domain_index, reverse=False)
			for v in val:
				domain_name_dict[v] = "%s|%s%s" % (v, domain_ab, domnr)
				domnr += 1
	elif strand_info_bgc_pos == 0:
		domain_name_dict = {}
		keys_sorted = sorted([int(k) for k in accessiondict.keys()], reverse = True)
		domnr = 1
		for key in keys_sorted:
			domain_index = [k.split('|')[-1].replace(domain_ab,'') for k in accessiondict.get(str(key))]
			val = sort_XonY(accessiondict.get(str(key)), domain_index, reverse=False)
			for v in val:
				domain_name_dict[v] = "%s|%s%s" % (v, domain_ab, domnr)
				domnr += 1
	elif strand_info_bgc_neg != 0 and strand_info_bgc_pos != 0:
		pattern1 = re.compile('\+\-')
		pattern2 = re.compile('\-\+')
		p1 = re.search(pattern1, strand_info)
		p2 = re.search(pattern2, strand_info)
		if (p1 == None and p2 != None) or (p1 != None and p2 == None):
			neg_strand = [k for k in accession if k.split('|')[3] == '-']
			pos_strand = [k for k in accession if k.split('|')[3] == '+']

			accession_common_neg = ["|".join(a.split("|")[0:3]) for a in neg_strand]
			accessiondict_neg = {}
			for x,y in zip(neg_strand, accession_common_neg):
				accessiondict_neg.setdefault(y.split('|')[2].split('-')[0],[]).append(x)

			domain_name_dict_neg = {}
			keys_sorted_neg = sorted([int(k) for k in accessiondict_neg.keys()], reverse = True)
			domnr_neg = 1
			for key in keys_sorted_neg:
				domain_index_neg = [k.split('|')[-1].replace(domain_ab,'') for k in accessiondict_neg.get(str(key))]
				val = sort_XonY(accessiondict_neg.get(str(key)), domain_index_neg, reverse=False)
				for v in val:
					if strand_info_bgc_neg >= strand_info_bgc_pos:
						domain_name_dict_neg[v] = "%s|%s%s" % (v, domain_ab, domnr_neg)
					else:
						domain_name_dict_neg[v] = "%s|%s%sb" % (v, domain_ab, domnr_neg)
					domnr_neg += 1

			accession_common_pos = ["|".join(a.split("|")[0:3]) for a in pos_strand]
			accessiondict_pos = {}
			for x,y in zip(pos_strand, accession_common_pos):
				accessiondict_pos.setdefault(y.split('|')[2].split('-')[0],[]).append(x)

			domain_name_dict_pos = {}
			keys_sorted_pos = sorted([int(k) for k in accessiondict_pos.keys()], reverse = False)
			domnr_pos = 1
			for key in keys_sorted_pos:
				domain_index_pos = [k.split('|')[-1].replace(domain_ab,'') for k in accessiondict_pos.get(str(key))]
				val = sort_XonY(accessiondict_pos.get(str(key)), domain_index_pos, reverse=False)
				for v in val:
					if strand_info_bgc_neg >= strand_info_bgc_pos:
						domain_name_dict_pos[v] = "%s|%s%sb" % (v, domain_ab, domnr_pos)
					else:
						domain_name_dict_pos[v] = "%s|%s%s" % (v, domain_ab, domnr_pos)
					domnr_pos += 1

			domain_name_dict = domain_name_dict_neg.copy()
			domain_name_dict.update(domain_name_dict_pos)

		else:
			domain_name_dict = {}
			for a in accession:
				domain_name_dict[a] = "%s|KS?" %(a)
	else:
		domain_name_dict = {}
		for a in accession:
			domain_name_dict[a] = "%s|KS?" %(a)


	return(domain_name_dict)


#--------------------------------
#-- run_index_domain
#--------------------------------
def run_index_domain(domain_id_list, domain_ab):
	strand_info = get_strand_info(domain_info = domain_id_list)
	new_domain_id = index_domain(accession=domain_id_list, strand_info=strand_info, domain_ab=domain_ab)
	return(new_domain_id)

