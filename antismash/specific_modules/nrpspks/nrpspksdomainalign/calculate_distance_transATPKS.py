import os
import re
import math
import pandas as pd
import itertools
import numpy as np
import logging
import munkres
import collections

from antismash.utils import sort_XonY
from datetime import datetime

def pair_domains_valid(A, B, BGCs):
    cluster_a = BGCs[A]
    cluster_b = BGCs[B]
    cluster_a_keys_1 = cluster_a.keys()
    cluster_b_keys_1 = cluster_b.keys()
    cluster_a_keys = []
    cluster_b_keys = []
    for ka in cluster_a_keys_1:
        if ka != 'clade_not_conserved' and ka != 'tree_not_conserved' and ka != 'out_group_like' and ka != 'noInfo':
            cluster_a_keys.append(ka)
    for kb in cluster_b_keys_1:
        if kb != 'clade_not_conserved' and kb != 'tree_not_conserved' and kb != 'out_group_like' and kb != 'noInfo':
            cluster_b_keys.append(kb)

    if(len(cluster_a_keys)>0 and len(cluster_b_keys)>0):
        return True
    else:
        return False

def readin_dist_matrix(dist_matrix):
    dist_matrix_input = [d for d in open(dist_matrix, 'r').read().split('\n') if d != '']
    index = [0, 1, 9]
    domain_name = []
    score_len = []
    for l in dist_matrix_input:
        s = len(l.split(','))
        score_len.append(s)   
    for l in dist_matrix_input:
        #fullheader = l.split(",")[0]
        #print fullheader
        #i0 = fullheader.split("|")[0]
        #i1 = fullheader.split("|")[1]
        #i9 = fullheader.split("|")[9]
        domain = "|".join([l.split(",")[0].split("|")[i] for i in index])
        domain_name.append(domain)
    domain_index = collections.defaultdict(list)
    for n, i in zip(domain_name, range(0, len(domain_name))):
        domain_index[n].append(i)
    return domain_index, dist_matrix_input

def __order_value_by_key(input_dict):
    out_list_1 = []
    out_list_2 = []
    if len(input_dict) != 0:
        for i in range(1, len(input_dict)+1):
            out_list_1.append(input_dict[i][0])
            out_list_2.append(input_dict[i][1])
    return out_list_1, out_list_2

def _get_new_cluster_index(ksa_per_new_cluster):
    new_cluster = []
    new_cluster_index = []
    for c in ksa_per_new_cluster.keys():
        cluster_i = ksa_per_new_cluster[c]
        cluster_id = '|'.join(cluster_i.keys()[0].split('|')[0:2])
        cluster_i_dict_wb = {}
        cluster_i_dict_wob = {}
        for k in cluster_i.keys():
            ksa_index = k.split('|')[9]
            index = re.sub('KSA', '', ksa_index)
            if 'b' in index:
                index_sub = re.sub('b', '', index)
                #cluster_i_dict_wb[int(index_sub)] = [ksa_index, cluster_i[k][0]]
                cluster_i_dict_wb[int(index_sub)] = [ksa_index, cluster_i[k]]
            else:
                #cluster_i_dict_wob[int(index)] = [ksa_index, cluster_i[k][0]]
                cluster_i_dict_wob[int(index)] = [ksa_index, cluster_i[k]]

        cluster_i_index_list_wb, cluster_i_list_wb = __order_value_by_key(input_dict = cluster_i_dict_wb)
        cluster_i_index_list_wob, cluster_i_list_wob = __order_value_by_key(input_dict=cluster_i_dict_wob)
        cluster_i_index = ','.join(cluster_i_index_list_wob+cluster_i_index_list_wb)
        cluster_i = ','.join(cluster_i_list_wob+cluster_i_list_wb)
        new_cluster.append('\n'.join([cluster_id, cluster_i]))
        new_cluster_index.append('\n'.join([cluster_id, cluster_i_index]))
    return new_cluster, new_cluster_index

def get_distance(compd1, compd2, index1, index2, score_lists, domain_index):
    """
    :param DM: DM is the pairwise distance matrix, where 0 represents two alignments are completely the same and 1 represents two alignment differ a lot
    :param compd1:
    :param compd2:
    :param index1:
    :param index2:
    :return: a list, where each entry is the sequence similarity score, 0 represents two alignments differs a lot and 1 represents two alignments are quite similar
    """
    A = ['%s|%s' % t for t in zip([compd1]*len(index1), index1)]
    B = ['%s|%s' % t for t in zip([compd2]*len(index2), index2)]
    output = []
    for i in range(len(A)):
        A_i = A[i]
        B_i = B[i]
        #row_index = domain_index.get(A_i)[0]
        #col_index = domain_index.get(B_i)[0]
        if A_i in domain_index:
            row_index = domain_index[A_i][0]
        else: #A not in domain_index
            output.append(-1)
            continue
        if B_i in domain_index:
            col_index = domain_index[B_i][0]
            #print "GOOD B"
        else: #B not in domain index
            output.append(-1)
            continue
        if row_index <= col_index:
            score_list = score_lists[col_index]
            col = row_index
        else:
            score_list = score_lists[row_index]
            col = col_index
        #score = round(1-float(score_list[col+1]),3)#similarity
        score = float(score_list[col+1])#frac_id
        output.append(score)
    return output

def generate_BGCs_DMS(cluster, cluster_KSindex, pseudo_aa, ksa_per_new_cluster, domain_index):
    new_cluster, new_cluster_index = _get_new_cluster_index(ksa_per_new_cluster=ksa_per_new_cluster)
    input = open(cluster, 'r').read().split('>')[1:]
    input = input + new_cluster
    cluster_id = [i.split('\n')[0] for i in input]
    pseudo_seq = [i.split('\n')[1] for i in input]
    cluster_seq = dict(zip(cluster_id, pseudo_seq))
    input_index = open(cluster_KSindex,'r').read().split('>')[1:]
    input_index = input_index + new_cluster_index
    index_cluster_id = [i.split('\n')[0] for i in input_index]
    index_pseudo_seq = [i.split('\n')[1] for i in input_index]
    cluster_seq_index = dict(zip(index_cluster_id, index_pseudo_seq))
    BGCs = {}
    """
    BGCs -- dictionary of this structure:  BGCs = {'cluster_name_x': { 'general_domain_name_x' : ['specific_domain_name_1', 'specific_domain_name_2'] } }
        - cluster_name_x: cluster name (can be anything)
        - general_domain_name_x: PFAM ID, for example 'PF00550'
        - specific_domain_name_x: ID of a specific domain that will allow to you to map it to names in DMS unequivocally
     (for example, 'PF00550_start_end', where start and end are genomic positions)."""
    for id in cluster_seq.keys():
        cluster_i = {}
        seq = cluster_seq[id].split(',')
        index = cluster_seq_index[id].split(',')
        seq_dict = {}
        for i in range(len(seq)):
            ks_id = index[i]
            seq_dict[id+';'+ks_id] = seq[i]
        domain_name = list(set(seq_dict.values()))
        for d in domain_name:
            key = d
            value = []
            for index, domain in seq_dict.items():
                if domain == d:
                    value.append(index)
            cluster_i[key] = value
        BGCs[id] = cluster_i
    DMS = {}
    """DMS -- dictionary of this structure: DMS = {'general_domain_name_x': { ('specific_domain_name_1',
    'specific_domain_name_2'): (sequence_identity, alignment_length), ... }   }
        - general_domain_name_x: as above
        - ('specific_domain_name_1', 'specific_domain_name_2'): pair of specific domains, sorted alphabetically
        - (sequence_identity, alignment_length): sequence identity and alignment length of the domain pair"""
    aa_list = [a for a in open(pseudo_aa, 'r').read().split('\n')[1:] if a != '']
    AA = [m.split('\t')[1] for m in aa_list]
    score_lists = [i.split(",") for i in DM_input]
    for aa in AA:
        domain_w_aa = []
        for i in BGCs.keys():
            cluster_i = BGCs[i]
            cluster_i_aa = cluster_i.keys()
            if aa in cluster_i_aa:
                domain_w_aa = domain_w_aa + cluster_i[aa]
        aa_pair_list = [tuple(sorted(list(i))) for i in list(itertools.combinations(domain_w_aa, 2))]
        aa_pair_dict = {}
        for p in aa_pair_list:
            index1 = [p[0].split(';')[1]]
            index2 = [p[1].split(';')[1]]
            compd1 = p[0].split(';')[0]
            compd2 = p[1].split(';')[0]
            score = get_distance(compd1, compd2, index1, index2, score_lists, domain_index)[0]
            aa_pair_dict[p] = (score, 0) ## Set alignment length to 0 always???
        DMS[aa] = aa_pair_dict
    return BGCs, DMS, cluster_seq, new_cluster, new_cluster_index

def calculate_GK(A, B, nbhood): #nbhood = 3, can be changed
    # calculate the Goodman-Kruskal gamma index
    #
    # example input:
    # A = ['IV_b', 'III_b', 'VI_c', 'II_g', 'II_c3', 'II_g', 'XV', 'I_a', 'VIII', 'III_b', 'clade_not_conserved']
    # B = ['clade_not_conserved', 'III_a', 'III_a', 'II_c3', 'III_a', 'III_a']
    GK = 0.
    if len(set(A) & set(B)) > 1:
        pairsA_1 = set([(A[i],j) for i in xrange(len(A)) for j in A[i+1:i+nbhood]])
        pairsB_1 = set([(B[i],j) for i in xrange(len(B)) for j in B[i+1:i+nbhood]])
        # exclude the domain pairs contain 'clade_not_conserved' or 'tree_not_conserved' or 'out_group_like' or 'noInfo'
        pairsA = []
        pairsB = []
        for a in pairsA_1:
            if 'clade_not_conserved' not in a and 'tree_not_conserved' not in a and 'out_group_like' not in a and 'noInfo' not in a:
                pairsA.append(a)
        for b in pairsB_1:
            if 'clade_not_conserved' not in b and 'tree_not_conserved' not in b and 'out_group_like' not in b and 'noInfo' not in b:
                pairsB.append(b)
        allPairs = set(pairsA + pairsB)
        Ns, Nr = 0.,0.
        for p in allPairs:
            if p in pairsA and p in pairsB: Ns += 1
            elif p in pairsA and tuple(p[::-1]) in pairsB: Nr += 1
            elif tuple(p[::-1]) in pairsA and p in pairsB: Nr += 1
            else: pass
        if (Nr + Ns) == 0:
            gamma = 0
        else:
            gamma = abs(Ns-Nr) / (Nr+Ns)
        GK = (1+gamma)/2.
    return GK

def cluster_distance(A, B, nbhood, BGCs, DMS, cluster_seq):
    # key is name of BGC, values is list of specific pfam domain names
    cluster_a = BGCs[A]  # will contain a dictionary where keys are pfam domains, and values are domains that map to a specific sequence in the DMS variable
    cluster_b = BGCs[B]  # { 'general_domain_name_x' : ['specific_domain_name_1', 'specific_domain_name_2'], etc }
    # A and B are lists of pfam domains
    try:
        # calculates the intersect
        cluster_a_keys_1 = cluster_a.keys()
        cluster_b_keys_1 = cluster_b.keys()
        cluster_a_keys = []
        cluster_b_keys = []
        for ka in cluster_a_keys_1:
            if ka != 'clade_not_conserved' and ka != 'tree_not_conserved' and ka != 'out_group_like' and ka != 'noInfo':
                cluster_a_keys.append(ka)
        for kb in cluster_b_keys_1:
            if kb != 'clade_not_conserved' and kb != 'tree_not_conserved' and kb != 'out_group_like' and kb != 'noInfo':
                cluster_b_keys.append(kb)

        jaccard = len(set(cluster_a_keys) & set(cluster_b_keys)) / \
                  float(len(set(cluster_a_keys)) + len(set(cluster_b_keys))
                        - len(set(cluster_a_keys) & set(cluster_b_keys)))
    except ZeroDivisionError:
        logging.error("Zerodivisionerror during the jaccard distance calculation. Can only happen when one or more clusters contains no pfam domains (or if those pfam domains are not from a conserved clade).")
        logging.debug("keys of cluster_a %s %s", A, cluster_a_keys)
        logging.debug("keys of cluster_b %s %s", B, cluster_b_keys)
        raise ValueError("One or more clusters contain no TransAT PFAM domains.")

    intersect = set(cluster_a.keys()).intersection(cluster_b.keys())  # shared pfam domains
    not_intersect = set(cluster_a.keys() + cluster_b.keys()) - intersect

    # diff_abundance: The difference in abundance of the domains per cluster
    # max_domains: Max occurence of each domain
    diff_abundance = 0
    max_domains = 0
    sum_distance = 0
    pair = ""
    for unshared_domain in not_intersect:
        # no need to look at seq identity or anchors, since these domains are unshared
        try:
            dom_set = cluster_a[unshared_domain]
        except KeyError:
            dom_set = cluster_b[unshared_domain]
        diff_abundance += len(dom_set)
        max_domains += len(dom_set)
    for shared_domain in intersect:
        seta = cluster_a[shared_domain]
        setb = cluster_b[shared_domain]
        if len(seta+setb) == 2:  # The domain occurs only once in both clusters
            pair = tuple(sorted([seta[0], setb[0]]))
            try:
                sum_distance = 1 - DMS[shared_domain][pair][0]
            except KeyError:
                print("KeyError on", pair)
                errorhandle = open(str(A)+".txt", 'w')
                print "\t"+str(pair)
                errorhandle.write(str(pair)+"\n")
                print "\t"+str(DMS[shared_domain])
                errorhandle.write(str(DMS[shared_domain]))
                errorhandle.close()
                print "\n"+str(DMS)
            max_domains += max(len(seta),len(setb))
            diff_abundance += sum_distance
        else:   # The domain occurs more than once in both clusters
            DistanceMatrix = [[1 for col in range(len(setb))] for row in range(len(seta))]
            for domsa in range(len(seta)):
                for domsb in range(len(setb)):
                    pair = tuple(sorted([seta[domsa], setb[domsb]]))
                    try:
                        Similarity = DMS[shared_domain][pair][0]
                    except KeyError:
                        print("KeyError on", pair)
                        errorhandle = open(str(B)+".txt", 'w')
                        errorhandle.write(str(pair)+"\n")
                        errorhandle.write(str(DMS[shared_domain]))
                        errorhandle.close()
                    seq_dist = 1 - Similarity
                    DistanceMatrix[domsa][domsb] = seq_dist
            # Only use the best scoring pairs
            Hungarian = munkres.Munkres()
            BestIndexes = Hungarian.compute(DistanceMatrix)
            accumulated_distance = sum([DistanceMatrix[bi[0]][bi[1]] for bi in BestIndexes])
            sum_distance = (abs(len(seta)-len(setb)) + accumulated_distance)  # diff in abundance + sequence distance
            max_domains += max(len(seta), len(setb))
            diff_abundance += sum_distance
    #  calculate the Goodman-Kruskal gamma index
    A_pseudo_seq = cluster_seq[A].split(',')
    B_pseudo_seq = cluster_seq[B].split(',')
    Ar = [item for item in A_pseudo_seq]
    Ar.reverse()
    GK = max([calculate_GK(A_pseudo_seq, B_pseudo_seq, nbhood = nbhood), calculate_GK(Ar, B_pseudo_seq, nbhood = nbhood)])
    diff_abundance /= float(max_domains)
    diff_abundance /= float(max_domains)
    diff_abundance = math.exp(-diff_abundance)  # transform from distance to similarity score
    Distance = 1 - (Jaccardw * jaccard) - (DDSw * diff_abundance) - (GKw * GK)
    Similarity_score = (Jaccardw * jaccard) + (DDSw * diff_abundance) + (GKw * GK)
    if Distance < 0:
        print("negative distance", Distance, "diff_abundance", diff_abundance, pair)
        print("Probably a rounding issue")
        print("Distance is set to 0 for these clusters")
        Distance = 0
    return Similarity_score, [A, B, jaccard, GK, diff_abundance]

def generate_distance_matrix(cluster_list, scale, nbhood, BGCs, DMS, cluster_seq):
    r = len(cluster_list)
    c = len(cluster_list)
    Dist = [[0 for x in range(r)] for y in range(c)]
    df = pd.DataFrame(Dist, index = cluster_list, columns = cluster_list)
    pairs = list(itertools.combinations(cluster_list, 2))
    detail_score = []
    for p in pairs:
        if(pair_domains_valid(p[1], p[0], BGCs)):
            sim_score, detail_score_p = cluster_distance(A = p[1],  B=p[0], nbhood=nbhood, BGCs=BGCs, DMS=DMS, cluster_seq=cluster_seq)
            rowname = p[1]
            colname = p[0]
            df.ix[rowname, colname] = 1 - sim_score/scale
            detail_score.append(detail_score_p)
    return df

# -------------------------------------------------------------------------------------------------
# -- function: average_dist
# -- This is the function to calculate the average distance between two groups of BGCs.
# -- In each group, the BGCs share the same composition of KS domains
# --------------------------------------------------------------------------------------------------
def average_dist(df_dist, unique_BGCs_info, len_assemblyline, ksa_per_new_cluster, BGCs):
    #-- please remove the BGCs in unique_BGCs_info file, i.e., KR857273|c1, KR857271|c1 and KR857272|c1, which are just replicated of misakinolid, luminaolid and tolytoxin
    new_cluster, _ = _get_new_cluster_index(ksa_per_new_cluster=ksa_per_new_cluster)
    new_bgc_id = [i.split('\n')[0] for i in new_cluster]
    colname = list(df_dist.dtypes.index)
    colname_index = dict(zip(colname, range(len(colname))))
    unique_BGCs = [u for u in open(unique_BGCs_info, 'r').read().split('\n') if u != '']
    unique_BGCs_list = []
    for b in unique_BGCs:
        v = b.split(': ')[1]
        k = b.split(': ')[0].split(',')
        if len(k) > len_assemblyline:
            unique_BGCs_list.append(v)
    unique_BGCs_list = unique_BGCs_list+[i.split('\n')[0] for i in new_cluster]
    # generate a list of new column names for the new distance matrix
    Dist_out_colname = []
    for u in unique_BGCs_list:
        u1 = u.split(',')
        if len(u1) == 1:
            Dist_out_colname.append(u)
        else:
            Dist_out_colname.append(u1[0]+'(extra '+str(len(u1)-1)+')')
    Dist_out_colname_dict = dict(zip(unique_BGCs_list, Dist_out_colname))
    Dist_out_colname_index_dict = dict(zip(Dist_out_colname, range(len(Dist_out_colname))))
    pairs_uniBGCs = list(itertools.combinations(unique_BGCs_list, 2))
    r = len(unique_BGCs_list)
    c = len(unique_BGCs_list)
    dist_out = [[0 for x in range(r)] for y in range(c)]
    Dist_out = pd.DataFrame(dist_out, index = Dist_out_colname, columns = Dist_out_colname)
    for p in pairs_uniBGCs:
        cluster1 = p[0].split(',')
        cluster2 = p[1].split(',')
        pairwise_cluster = list(itertools.product(cluster1, cluster2))
        dist_per_pair = []
        for pair in pairwise_cluster:
            if(pair[1] in colname_index.keys() and pair[0] in colname_index.keys()):
                col_index = colname_index[pair[1]]
                row_index = colname_index[pair[0]]
                if col_index < row_index:
                    dist_per_pair.append(df_dist.loc[pair[0], pair[1]])
                else:
                    dist_per_pair.append(df_dist.loc[pair[1], pair[0]])
        if(len(dist_per_pair)>0):
            ave = round(np.mean(dist_per_pair), 3)
        else:
            ave=0 # This needs to be double checked...Assuming now that 0 means very dissimilar
        ave_rowname = Dist_out_colname_dict[p[0]]
        ave_colname = Dist_out_colname_dict[p[1]]
        ave_rowname_index = Dist_out_colname_index_dict[ave_rowname]
        ave_colname_index = Dist_out_colname_index_dict[ave_colname]
        if ave_rowname_index > ave_colname_index:
            Dist_out.loc[ave_rowname, ave_colname] = ave
        else:
            Dist_out.loc[ave_colname, ave_rowname] = ave
    return Dist_out, Dist_out_colname, new_bgc_id

def get_most_similar_bgcs(bgc_id_list, ave_dist_df, ave_dist_colname, cut_off_nr):
    r = ave_dist_df.shape[0]
    c = ave_dist_df.shape[1]
    for i in range(0,r):
        for j in range(0,c):
            if i == j:
                ave_dist_df.iloc[i,j] = 1
            elif i < j:
                ave_dist_df.iloc[i,j] = ave_dist_df.iloc[j,i]
    output_dict = {}
    for id in bgc_id_list:
        dist_score = list(ave_dist_df.loc[id, :])
        sorted_bgc = sort_XonY(ave_dist_colname, dist_score, reverse=False)
        top_ranked_bgc = sorted_bgc[0:cut_off_nr]
        top_ranked_bgc.append(id)
        pairs = list(itertools.permutations(top_ranked_bgc, 2))
        r_bgc = len(top_ranked_bgc)
        c_bgc = len(top_ranked_bgc)
        dist_out_bgc = [[0 for x in range(r_bgc)] for y in range(c_bgc)]
        Dist_out_bgc = pd.DataFrame(dist_out_bgc, index=top_ranked_bgc, columns=top_ranked_bgc)
        for p in pairs:
            Dist_out_bgc.loc[p[0], p[1]] = ave_dist_df.loc[p[0], p[1]]
        for i in range(0, r_bgc):
            for j in range(0, c_bgc):
                if i == j:
                    Dist_out_bgc.iloc[i,j] = 0.0
        output_dict[id] = [top_ranked_bgc, Dist_out_bgc]
    return output_dict

def run_calculate_distance(data_dir, seq_simScore, ksa_per_new_cluster, cutoff_bgc_nr):
    logging.info("TRANSATPKS: generating distance matrix of assembly lines")
    global Jaccardw, GKw, DDSw
    global domain_index, DM_input
    ## Set weights
    Jaccardw = 0.5
    GKw = 0.25
    DDSw = 0.25
    ## Set paths
    cluster = os.path.join(data_dir, 'completeBGCs_mc.faa')
    cluster_KSindex = os.path.join(data_dir, 'completeBGCs_mc_index.faa')
    pseudo_aa = os.path.join(data_dir, 'transAT_mc.txt')
    unique_BGCs_info = os.path.join(data_dir, "unique_mc.txt")
    ## Calc
    len_assemblyline = 2
    domain_index, DM_input = readin_dist_matrix(dist_matrix=seq_simScore)
    BGCs, DMS, cluster_seq, new_cluster, new_cluster_index = generate_BGCs_DMS(cluster, cluster_KSindex, pseudo_aa, ksa_per_new_cluster, domain_index)
    ## Structure of BGCs
    ##     dCluster_name -> dClade_name -> list_of_matching_KSs (format = cluster;KS)
    ##
    ## Structure of DMS
    ##     dClade_name -> d(pair -> (seqid, alen)
    ##
    #print "BGCs:\t"+str(BGCs)+"\n"
    #print "DMS:\t"+str(DMS)+"\n"
    #print "cluster_seq:\t"+str(cluster_seq)+"\n"
    #print "new_cluster:\t"+str(new_cluster)+"\n"
    #print "new_cluster_index:\t"+str(new_cluster_index)+"\n"
    print str(datetime.now())+"\tgenerate_distance_matrix"
    dist_score_assembly_line = generate_distance_matrix(cluster_list=cluster_seq.keys(), scale=1, nbhood=3, BGCs=BGCs, DMS=DMS, cluster_seq=cluster_seq) # currently pretty slow, potential speedup could be limiting this to assembly lines that have a set # of domains with clade overlap and/or moving to cython
    print str(datetime.now())+"\taverage_dist"
    ave_dist, dist_out_colname, new_bgc_id = average_dist(df_dist = dist_score_assembly_line, unique_BGCs_info=unique_BGCs_info, len_assemblyline=len_assemblyline, ksa_per_new_cluster=ksa_per_new_cluster, BGCs=BGCs) # this is also pretty slow
    print str(datetime.now())+"\tget_most_similar_bgcs"
    most_similar_bgcs = get_most_similar_bgcs(bgc_id_list=new_bgc_id, ave_dist_df=ave_dist, ave_dist_colname=dist_out_colname, cut_off_nr=cutoff_bgc_nr)
    return most_similar_bgcs, new_cluster, new_cluster_index
