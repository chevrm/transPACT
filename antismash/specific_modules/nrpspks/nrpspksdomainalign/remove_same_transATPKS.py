#!/usr/bin/env python
#-------------------------------------------------------------------------------------------
## Copyright (c) 2016 Xiaowen Lu
## Wageningen University
## Bioinformatics Group
#
#  This is the function to remove the assembly lines that have the same domain composition
#-------------------------------------------------------------------------------------------


import re
import argparse

'''
training_list = ['dorrigocin_migrastatin|GQ274953']
novel_list = ['CXWA01000005|c1','CGFR01000002|c1', 'JZDI01000001|c2', 'AWQY01000002|c2', 'AJST01000001|c6']
ID_list = ['AMPK01000010|c2', 'bacillaene_Bamy|AJ634060']
outfile_repeatID = '/Users/Xiaowen/Desktop/repeatID.txt'
annotateKS_training = '/Users/Xiaowen/Documents/trans_AT/trans_AT_LuX_version3/annotateKS_perCompound_v3_uniPK.txt'
annotateKS_novel = '/Users/Xiaowen/Documents/trans_AT/KS_domains_PredtransAT_output_20160330/annotateKS_per_pred_transATPKS.txt'
'''

#-------------------------------------
#--Function: extract_KS_anno_unique
#-------------------------------------
def extract_KS_anno_unique(ID_list, annotateKS_training, annotateKS_novel, outfile_repeatID):
	'''

	:param ID_list:
	:param annotateKS_training:
	:param annotateKS_novel:
	:return: a dict where key is

	'''

	annoKS_t = [t for t in open(annotateKS_training, 'r').read().split('\n') if t != '']
	annoKS_n1 = [t for t in open(annotateKS_novel, 'r').read().split('\n') if t != '']
	annoKS_n2 = [re.sub('clade_not_conserved', 'noInfo', n) for n in annoKS_n1]
	annoKS_n = [re.sub('out_group_like', 'noInfo', n) for n in annoKS_n2]
	annoKS = annoKS_t[1:] + annoKS_n[1:]
	KS_id = annoKS_n[0]

	output_anno = {}
	output_index = {}

	for i in ID_list:

		i_adj = i  + '\t'
		compd_i = [a for a in annoKS if i_adj in a][0]
		compd_id = compd_i.split('\t')[0]

		compd_ks_index = zip(compd_i.split('\t')[1:], KS_id.split('\t')[1:])
		compd_ks = [x for x,y in compd_ks_index if x != 'NA' ][1:]
		compd_index = [y for x,y in compd_ks_index if x != 'NA'][1:]

		output_anno[compd_id] = ';'.join(compd_ks)
		output_index[compd_id] = compd_index


	#-- the assembly lines are filtered only based on the domian composition, the domain sequence similarity are not considered
	output_anno_val = list(set(output_anno.values()))
	output_anno_reverse = {}

	for anno in output_anno_val:
		ID = []
		for x,y in output_anno.items():
			if y == anno:
				ID.append(x)

		output_anno_reverse[anno] = ID

	output_anno_mergeID = {}
	repeatID_dict = {}
	output_index_mergeID = {}
	for k in output_anno_reverse.keys():
		val = output_anno_reverse[k]
		if len(val) == 1:
			new_key = val[0]
			output_anno_mergeID[new_key] = k.split(';')
			output_index_mergeID[new_key] = output_index[new_key]
		else:
			new_key = ''.join([val[0], '(extra ', str(len(val) -1), ')'])
			repeatID_dict[new_key] = ';'.join(val)
			output_anno_mergeID[new_key] = k.split(';')
			output_index_mergeID[new_key] = output_index[val[0]]


	out = open(outfile_repeatID, 'w')
	for x,y in repeatID_dict.items():
		out.write('%s\t%s\n' %(x,y))
	out.close()

	return output_anno_mergeID, output_index_mergeID



# KS_anno_mergeID, KS_index_mergeID = extract_KS_anno_unique(ID_list = ID_list, annotateKS_training = annotateKS_training, annotateKS_novel = annotateKS_novel, outfile_repeatID = outfile_repeatID)

if __name__ == "__main__":
	"Run this file from here as a script"
	#Check if parameters are provided; if not, exit with explanation
	parser = argparse.ArgumentParser(add_help = True)
	parser.add_argument('--ID_list', nargs = '*')
	parser.add_argument('--annotateKS_training')#the input here can be a motif or a list of BGC id
	parser.add_argument('--annotateKS_novel')
	parser.add_argument('--outfile_repeatID')

	args = parser.parse_args()
	if any(t == [] for t in vars(args).values()) == False:
		ID_list = list(args.ID_list)
		annotateKS_training = args.annotateKS_training
		annotateKS_novel = args.annotateKS_novel
		outfile_repeatID = args.outfile_repeatID

		KS_anno_mergeID, KS_index_mergeID = extract_KS_anno_unique(ID_list = ID_list, annotateKS_training = annotateKS_training, annotateKS_novel = annotateKS_novel, outfile_repeatID = outfile_repeatID)

		for k in KS_anno_mergeID.keys():
			print('%s: %s\n' %(k, KS_anno_mergeID[k]))
		for k in KS_index_mergeID.keys():
			print('%s: %s\n' %(k, KS_index_mergeID[k]))

	else:
		print("Please provide needed profiles")


#python /Users/Xiaowen/Documents/trans_AT/scripts/remove_same_transATPKS.py --ID_list 'dorrigocin_migrastatin|GQ274953' 'CXWA01000005|c1' 'CGFR01000002|c1' 'JZDI01000001|c2' 'AWQY01000002|c2' 'AJST01000001|c6' --annotateKS_training /Users/Xiaowen/Documents/trans_AT/trans_AT_LuX_version3/annotateKS_perCompound_v3_uniPK.txt --annotateKS_novel /Users/Xiaowen/Documents/trans_AT/KS_domains_PredtransAT_output_20160330/annotateKS_per_pred_transATPKS.txt --outfile_repeatID /Users/Xiaowen/Desktop/repeatID.txt > /Users/Xiaowen/Destkop/remove_same_transATPKS.out
