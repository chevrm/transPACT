## Copyright (c) 2016 Xiaowen Lu
## Wageningen University
## Bioinformatics Group
##
##
#  This module is used to get two pieces of information,
#  i)get info. where the KS domains are located on the same gene
#  ii)whether the last KS is located at the end of the assembled genome, especially for the cases where only draft genome is available



# --import libs and add path to sys path
import re


# ---------------------------------------------------------------------------
# - function: get_module_loc_perBGC
# -
# - This function to get information on whether multiple KS domains are located
# - on the same gene
# ---------------------------------------------------------------------------
def get_module_loc_perBGC(KS_info, BGC_id, feature_loc, index):
    '''
    :param KS_info: read in '/Users/Xiaowen/Documents/trans_AT/KS_domains_PredtransAT_output_20160330/KSinCompound_CladeInfo_pred_transATPKS.txt' as list, where each line is a entry of the list
    :param BGC_id: e.g. CP000076|c4
    :param feature_loc: a numeric number, e.g., 6, which indicating the position of 'KND26274' in JPPY01000213|c1|21971-41386|-|IQ63_37575|polyketide_synthase|KND26274|7-435|KS1|KS1, counting from 0
    :param index: a list of ordered KS domains, e.g., ['KS1', 'KS2', 'KS3', 'KS4', 'KS1b', 'KS2b']
    :return: a list, ['F', 'T', 'F', 'F', 'F'], 'T' represents that a small black line needs to add between the domain and the domain next to it
    '''
    bgc_id_1 = ''.join(['\t', BGC_id.split('|')[0], '\t'])
    bgc_id_2 = ''.join(['\t', BGC_id.split('|')[1], '\t'])
    KS_list = [k for k in KS_info if bgc_id_1 in k and bgc_id_2 in k]

    KS_dict = {}  # this dict has a format as {'KS4': 'EEP84582', 'KS1b': 'EEP84596', 'KS1': 'EEP84584', 'KS3': 'EEP84583', 'KS2': 'EEP84583'}

    for ks in KS_list:
        v = ks.split('\t')[1].split('|')[feature_loc]
        k = ks.split('\t')[3]
        KS_dict[k] = v

    KS_loc_list = []  # put the gene info of each KS domain in the same order as the assembly: ['EEP84584', 'EEP84583', 'EEP84583', 'EEP84582', 'EEP84596']
    for i in index:
        KS_loc_list.append(KS_dict[i])

    output_boolean = []  # ['F', 'T', 'F', 'F', 'F'], 'T' represents that a small black line needs to add between the domain and the domain next to it
    for i in range(len(KS_loc_list) - 1):
        if KS_loc_list[i] == KS_loc_list[i + 1]:
            output_boolean.append('T')
        else:
            output_boolean.append('F')

    output_boolean.append('F')

    return output_boolean


# -------------------------------------------------------------------------------
# - function: get_contig_edge_perBGC
# -
# - This function to get information on the KS ids that are locted on the left/right
# - boarder of the assembled genome
# --------------------------------------------------------------------------------
def get_contig_edge_perBGC(KS_info, BGC_id, index, edge_len, feature_edge, feature_strand, genome_size):
    '''
    this function is needed only when an assembly line that is identified on a draft genome
    :param KS_info: read in '/Users/Xiaowen/Documents/trans_AT/KS_domains_PredtransAT_output_20160330/KSinCompound_CladeInfo_pred_transATPKS.txt' as list, where each line is a entry of the list
    :param BGC_id: e.g. CP000076|c4
    :param index: a list of ordered KS domains, e.g., ['KS1', 'KS2', 'KS3', 'KS4', 'KS1b', 'KS2b']
    :param edge_len: a fraction, such as 0.01, and edge_border = egde_len * genome_size
    :param feature_edge: a numeric number, e.g., 2, which indicating the position of '21971-41386' in JPPY01000213|c1|21971-41386|-|IQ63_37575|polyketide_synthase|KND26274|7-435|KS1|KS1, counting from 0
    :param feature_strand: a numeric number, e.g., 3, which indicating the position of '-' in JPPY01000213|c1|21971-41386|-|IQ63_37575|polyketide_synthase|KND26274|7-435|KS1|KS1, counting from 0
    :param genome_size: a numeric number of the size of the assembled sequence in a genbank file
    :return: two list, [KS1] contains KS id that is on the left border of the assembled genome, [KS5] contains KS id that is on the right border of the assembled genome
    '''

    left_bd = [1, 1 + edge_len]
    right_bd = [genome_size - edge_len, genome_size]

    # get the KS on the border of the assembly line
    index_lists1 = [[i for i in index if 'b' not in i], [i for i in index if 'b' in i]]
    index_lists = [l for l in index_lists1 if l != []]

    border_ks = []
    for l in index_lists:
        border_ks.append(l[0])
        border_ks.append(l[-1])

    # get genome locus_tag for each border KS
    bgc_id_1 = ''.join(['\t', BGC_id.split('|')[0], '\t'])
    bgc_id_2 = ''.join(['\t', BGC_id.split('|')[1], '\t'])
    KS_list = [k for k in KS_info if bgc_id_1 in k and bgc_id_2 in k]

    KS_dict = {}  # this dict has a format as {'KS4': 'EEP84582', 'KS1b': 'EEP84596', 'KS1': 'EEP84584', 'KS3': 'EEP84583', 'KS2': 'EEP84583'}
    KS_strand_dict = {}
    for ks in KS_list:
        v1 = ks.split('\t')[1].split('|')[feature_edge]
        v2 = ks.split('\t')[1].split('|')[feature_strand]
        k = ks.split('\t')[3]
        KS_dict[k] = v1
        KS_strand_dict[k] = v2

    border_ks_loc = {}
    for k in border_ks:
        border = KS_dict[k]
        border_ks_loc[k] = [int(b) for b in border.split('-')]

    border_info = {}
    for k in border_ks_loc.keys():
        start_loc_k = border_ks_loc[k][0]
        end_loc_k = border_ks_loc[k][1]
        strand_k = KS_strand_dict[k]

        v = ''
        if start_loc_k >= left_bd[0] and start_loc_k <= left_bd[1]:
            if strand_k == '-':
                v = 'R'
            else:
                v = 'L'
        elif end_loc_k >= right_bd[0] and end_loc_k <= right_bd[1]:
            if strand_k == '-':
                v = 'L'
            else:
                v = 'R'

        border_info[k] = v

    border_info_w_value = {}
    for k in border_info.keys():
        if border_info[k] != '':
            border_info_w_value[k] = border_info[k]

    KS_wo_b = [i for i in border_info_w_value.keys() if 'b' not in i]
    KS_w_b = [i for i in border_info_w_value.keys() if 'b' in i]

    left_break = []
    right_break = []
    if len(KS_wo_b) == 1:
        if border_info[KS_wo_b[0]] == 'R':
            right_break.append(KS_wo_b[0])
        elif border_info[KS_wo_b[0]] == 'L':
            left_break.append(KS_wo_b[0])

    elif len(KS_wo_b) == 2:
        KS_index = [int(re.sub('KS(A|)', '', k)) for k in KS_wo_b]
        KS_sorted = [x for (y, x) in sorted(zip(KS_index, KS_wo_b))]

        if border_info[KS_wo_b[0]] == 'R' and border_info[KS_wo_b[1]] == 'R':
            right_break.append(KS_sorted[1])
        elif border_info[KS_wo_b[0]] == 'L' and border_info[KS_wo_b[1]] == 'L':
            left_break.append(KS_sorted[0])
        elif border_info[KS_wo_b[0]] == 'R' and border_info[KS_wo_b[1]] == 'L':
            right_break.append(KS_wo_b[0])
            left_break.append(KS_wo_b[1])
        elif border_info[KS_wo_b[0]] == 'L' and border_info[KS_wo_b[1]] == 'R':
            right_break.append(KS_wo_b[1])
            left_break.append(KS_wo_b[0])

    if len(KS_w_b) == 1:
        if border_info[KS_w_b[0]] == 'R':
            right_break.append(KS_w_b[0])
        elif border_info[KS_w_b[0]] == 'L':
            left_break.append(KS_wo_b[0])
    elif len(KS_w_b) == 2:
        KS_w_index1 = [re.sub('KS(A|)', '', k) for k in KS_w_b]
        KS_w_index = [int(re.sub('b', '', k)) for k in KS_w_index1]
        KS_w_sorted = [x for (y, x) in sorted(zip(KS_w_index, KS_w_b))]

        if border_info[KS_w_b[0]] == 'R' and border_info[KS_w_b[1]] == 'R':
            right_break.append(KS_w_sorted[1])
        elif border_info[KS_w_b[0]] == 'L' and border_info[KS_w_b[1]] == 'L':
            left_break.append(KS_w_sorted[0])
        elif border_info[KS_w_b[0]] == 'R' and border_info[KS_w_b[1]] == 'L':
            right_break.append(KS_w_b[0])
            left_break.append(KS_w_b[1])
        elif border_info[KS_w_b[0]] == 'L' and border_info[KS_w_b[1]] == 'R':
            right_break.append(KS_w_b[1])
            left_break.append(KS_w_b[0])

    return left_break, right_break


# ---------------------------------------
# - function: run_get_module_loc_perBGC
# ---------------------------------------
def run_get_moduleLoc_contigEdge(KSindex, KS_info_file, feature_loc, genomeSize_file, feature_edge, feature_strand,
                                 edge_len, ksa_per_new_cluster, genomeSize_new_bgc, new_bgc_id):
    '''
    :param KSindex: a dict, {'CP000076|c4(extra 1)': ['KS1', 'KS2', 'KS3], 'bgc_id_2': ['KS1', 'KS2'],...}
    :param KS_info_file: a list containing the filenames, one input example, ['/Users/Xiaowen/Documents/trans_AT/KS_domains_PredtransAT_output_20160330/KSinCompound_CladeInfo_pred_transATPKS.txt']
    :param feature_loc: a numeric number, e.g., 6, which indicating the position of 'KND26274' in JPPY01000213|c1|21971-41386|-|IQ63_37575|polyketide_synthase|KND26274|7-435|KS1|KS1, counting from 0
    :param genomeSize_file: a list containing the filenames, e.g., ['/Users/Xiaowen/Documents/trans_AT/KS_domains_PredtransAT_output_20160330/genome_size_pred_transATPKS.txt']
    :param edge_len: a number (bp) indicating the border length
    :param feature_edge: a numeric number, e.g., 2, which indicating the position of '21971-41386' in JPPY01000213|c1|21971-41386|-|IQ63_37575|polyketide_synthase|KND26274|7-435|KS1|KS1, counting from 0
    :param feature_strand: a numeric number, e.g., 3, which indicating the position of '-' in JPPY01000213|c1|21971-41386|-|IQ63_37575|polyketide_synthase|KND26274|7-435|KS1|KS1, counting from 0
    :param bgc_id_pos: a numeric number to indication the position of the genbank id in a bgc_id, 'CP000076|c4'. In this case, bgc_id_pos = 0
    :return: a dict, {'CP000076|c4(extra 1)': [['T', 'T', 'T', 'F', 'T', 'T', 'F', 'T', 'T', 'F', 'T', 'T', 'F', 'T', 'F', 'F'], [KS id on the left border of the assembled genome], [KS id on the left border of the assembled genome]], ...}
    '''
    new_cluster_ks_info = []
    start_index = 0
    for k in ksa_per_new_cluster.keys():
        cluster_k = ksa_per_new_cluster[k]
        for ck in cluster_k.keys():
            v1 = ck
            v2 = ck.split('|')[0]
            v3 = ck.split('|')[9]
            v4 = ck.split('|')[1]
            v5 = cluster_k[ck][0]
            new_cluster_ks_info.append('\t'.join((str(start_index), v1, v2, v3, v4, v5)))
            start_index += 1

    ID = KSindex.keys()
    KS_info = []
    for f in KS_info_file:
        KS_info_p = [l for l in open(f, 'r').read().split('\n')[1:] if l != '']
        KS_info = KS_info + KS_info_p
    KS_info = KS_info + new_cluster_ks_info

    genomeSize_info = []
    for f in genomeSize_file:
        genomeSize_info_p = [g for g in open(f, 'r').read().split('\n') if g != '']
        genomeSize_info = genomeSize_info + genomeSize_info_p

    genomeSize_dict = {}
    for gs in genomeSize_info:
        k = gs.split('\t')[0]
        v = gs.split('\t')[1]
        genomeSize_dict[k] = v
    for new_id in new_bgc_id:
        genomeSize_dict[new_id] = str(genomeSize_new_bgc)

    output = {}
    for id in ID:
        id_1 = id.split('(')[0]
        if len(id_1.split('|')[1]) > 4:
            bgc = id_1.split('|')[1]
        else:
            bgc = id_1.split('|')[0]

        index_id = KSindex[id]
        moduleLoc_id = get_module_loc_perBGC(KS_info=KS_info, BGC_id=id_1, feature_loc=feature_loc, index=index_id)

        genome_size = int(genomeSize_dict[id_1])
        if len(bgc) == 8:
            KSID_left = []
            KSID_right = []
        else:
            KSID_left, KSID_right = get_contig_edge_perBGC(KS_info=KS_info, BGC_id=id_1, index=index_id,
                                                           feature_edge=feature_edge, genome_size=genome_size,
                                                           feature_strand=feature_strand, edge_len=edge_len)

        output[id] = [moduleLoc_id, KSID_left, KSID_right]

    return output


'''
KSindex = KS_index
KSindex = KS_index_mergeID
edge_len = 10000 #in bp, indicate the length of border on a contig
genomeSize_file = '/Users/Xiaowen/Documents/trans_AT/KS_domains_PredtransAT_output_20160330/genome_size_pred_transATPKS.txt'
KS_info_file = '/Users/Xiaowen/Documents/trans_AT/KS_domains_PredtransAT_output_20160330/KSinCompound_CladeInfo_pred_transATPKS.txt'

genomeSize_file = '/Users/Xiaowen/Documents/trans_AT/KS_domains_PredtransAT_output_20160330/genome_size_training_transATPKS.txt'
KS_info_file = '/Users/Xiaowen/Documents/trans_AT/trans_AT_LuX_version3/KSinCompound_CladeInfo_v4_exclude.txt'
feature_loc = 6
feature_edge = 2
feature_strand = 3


domain_loc_edge = run_get_moduleLoc_contigEdge(KSindex = KSindex, KS_info_file = KS_info_file, feature_loc = feature_loc, genomeSize_file = genomeSize_file, feature_edge = feature_edge, feature_strand = feature_strand, edge_len = edge_len)
'''
