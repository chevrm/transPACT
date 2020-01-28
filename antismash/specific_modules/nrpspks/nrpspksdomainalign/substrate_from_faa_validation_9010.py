import logging
import os
import re
import subprocess
import copy
from ete2 import Tree
from itertools import islice
from Bio.Align.Applications import MuscleCommandline

from Bio import SeqIO
import sys
sys.path.append(os.path.abspath('../../../..')) ## just for testing

from antismash import utils
#from antismash.utils import sort_XonY
import time

#UTILITY
def reformat(input, out_filename):
    """
    :Description: This is the function to transform the aligned fasta file from multiple line version to single line version, which is usually the input format of mega

    :param input: read in a fasta file
    :param out_filename: a string, which is name of the output file
    :return: NA. the output is saved as out_filename
    """

    entries = input.split('>')[1:]
    out = open(out_filename, 'w')

    for l in entries:
        id = re.sub(r'(\:|\'|\(|\)|\,|\?|\;)', '', l.split('\n')[0])
        seq = ''.join(l.split('\n')[1:])
        out.write('>%s\n%s\n' % (id, seq))

    out.close()


###REFACTORED PART BELOW
def align_ks_domains(reference_alignment, ks_names, ks_seqs, data_dir):
    """Function that aligns a number of query KS domain sequences to the 
    reference alignment of KS domains.
    """
    #Set file names and write query domains to temp input file
    in_temp = os.path.join(os.getcwd(), "in_seq.fasta")
    in_temp_aligned = os.path.join(os.getcwd(), "in_seq_aligned.fasta")
    out_temp = os.path.join(os.getcwd(), "out_seq.fasta")
    alignment_file = os.path.join(os.getcwd(), "aligned.fasta")
    with open(in_temp, "w") as tmp_input:
        for name, seq in zip(ks_names, ks_seqs):
            tmp_input.write("%s\n%s\n" % (name, seq))

    #Generate alignment of query sequences
    muscle_cmd = str(MuscleCommandline(input=in_temp, out=in_temp_aligned))
    out, err, retcode = utils.execute(muscle_cmd.split(" "))
    if retcode == 1:
        logging.error("Alignment of query KS sequences with Muscle failed. Check if Muscle is installed appropriately.")
        sys.exit(1)

    #Align the query alignment to the reference alignment using muscle --profile
    muscle_cmd = str(MuscleCommandline(profile='True', in1=reference_alignment, in2=in_temp_aligned, out=out_temp))
    out, err, retcode = utils.execute(muscle_cmd.split(" "))
    if retcode == 1:
        logging.error("Alignment of query+reference KS sequences with Muscle failed. Check if Muscle is installed appropriately.")
        sys.exit(1)
    else:
        f_temp_input = open(out_temp, 'r').read()
        reformat(input=f_temp_input, out_filename=alignment_file)

    #Remove temporary files
    for f in [in_temp, out_temp]:
        os.remove(f)

    return alignment_file


def run_pplacer(reference_alignment, alignment_file, data_dir, reference_tree):
    """Function that uses the reference tree with the new alignment to place
    query domains onto reference tree.
    """
    #Locations of files
    #reference_tree = os.path.join(data_dir, "RAxML_bestTree.647KS_RAxML.tre")
    pplacer_json = os.path.join(os.getcwd(), "pplacer_tree.jplace")

    #Reference package creation: taxit create --aln-fasta test_set_for_development.fasta --tree-stats test_set_for_development.log --tree-file test_set_for_development.nwk -P test_set_for_development.refpkg -l test
    #Reference package creation: taxit create --aln-fasta 647KS_mcformat.afa --tree-stats RAxML_info.647KS_RAxML.tre --tree-file RAxML_bestTree.647KS_RAxML.tre -P RAxML_bestTree.647KS_RAxML.refpkg -l 647KS
    
    pplacer_cmd = ["pplacer", "-t", reference_tree, "-r", reference_alignment, "-o", pplacer_json, "-c", os.path.join(data_dir, "RAxML_bestTree.647KS_RAxML.refpkg"), alignment_file]
    out, err, retcode = utils.execute(pplacer_cmd)
    if retcode == 1:
        logging.error("Running pplacer failed. Check if the program is installed appropriately.")
        sys.exit(1)
    guppy_cmd = ["guppy", "sing", pplacer_json]
    out, err, retcode = utils.execute(guppy_cmd)
    if retcode == 1:
        logging.error("Running guppy failed. Check if the program is installed appropriately.")
        sys.exit(1)
    return os.getcwd() + os.sep + "pplacer_tree.sing.tre"

def get_leaf2clade(leaf2cladetbl):
    leaf2clade = {}
    clade2ann = {}
    with open(leaf2cladetbl) as c:
        for ln in c.read().splitlines():
            ksname, clade, ann = ln.split("\t")
            leaf2clade[ksname] = clade
            clade2ann[clade] = ann
    return(leaf2clade, clade2ann)

def get_clade(q, tree_hits, t, funClades):
    cladelist = {}
    qname = tree_hits[q].name
    newtree = copy.deepcopy(t)
    pref = re.sub(r"^(.+)_#\d+_M=\d+?\.?\d*$", "\g<1>", qname)
    keep = []
    lcount = 0
    for leaf in newtree.get_leaves():
        lcount += 1
        ln = leaf.name.split("_")
        if re.match("^#\d+$", ln[-2]) is None:
            keep.append(leaf)
        elif ln[-2] == '#'+q:
            keep.append(leaf)
    newtree.prune(keep)
    grandparent = None
    for leaf in newtree.get_leaves():
        if leaf.name == qname:
            grandparent = leaf.up.up
            break
    for leaf in grandparent.get_leaves():
        if leaf.name == qname:
            continue
        elif leaf.name in funClades:
            if funClades[leaf.name] == 'query_seq':
                continue
            elif funClades[leaf.name] in cladelist.keys():
                cladelist[funClades[leaf.name]] += 1
            else:
                cladelist[funClades[leaf.name]] = 1
        else:
            if 'unknown_clade' in cladelist.keys():
                cladelist['unknown_clade'] += 1
            else:
                cladelist['unknown_clade'] = 1
    clade_assignment = 'clade_not_conserved'
    if len(cladelist.keys()) == 1:
        for k in cladelist:
            clade_assignment = k
    return pref, clade_assignment

def parse_pplacer(pplacer_tree, data_dir, masscutoff, mode):
    #leaf2cladetbl = data_dir + '/clade_mc.csv'
    leaf2cladetbl = data_dir + '/annotation_mc_20180727.tsv'
    ## Get clade designations
    funClades, clade2ann = get_leaf2clade(leaf2cladetbl)
    monoclade = {}
    with open(pplacer_tree) as f:
        for ln in f.read().splitlines():
            t = Tree(ln)
            tree_hits = {}
            ## Identify the list of pplacer placements
            for leaf in t.get_leaves():
                ln = leaf.name.split("_")
                if re.match("^#\d+$", ln[-2]) is not None:
                    n = re.sub(r"^#(\d+)$", "\g<1>", ln[-2])
                    tree_hits[n] = leaf
                    funClades[leaf.name] = 'query_seq'
            ## Look to see when threshold is met
            if mode == 'top_down':
                to_check = []
                running_mass = 0
                for placement in sorted(tree_hits): ## NOTE: unsure about behavior when n>=10
                    running_mass += float(re.sub(r"^.+#\d+_M=(\d+?\.?\d*)$", "\g<1>", tree_hits[placement].name))
                    to_check.append(placement)
                    if(running_mass >= masscutoff):
                        break
                cl = {}
                for q in to_check:
                    pref, clade_assignment = get_clade(q, tree_hits, t, funClades)
                    if clade_assignment in cl.keys():
                        cl[clade_assignment] += 1
                    else:
                        cl[clade_assignment] = 1;
                if len(cl.keys()) == 1:
                    for k in cl:
                        monoclade[pref] = k
                else:
                    monoclade[pref] = 'clade_not_conserved'
            elif mode == 'additive':
                totalmass = {}
                pref = ''
                for placement in tree_hits:
                    mass = float(re.sub(r"^.+#\d+_M=(\d+?\.?\d*)$", "\g<1>", tree_hits[placement].name))
                    pref, clade_assignment = get_clade(placement, tree_hits, t, funClades)
                    if clade_assignment in totalmass:
                        totalmass[clade_assignment] += mass
                    else:
                        totalmass[clade_assignment] = mass
                best_clade = max(totalmass, key=totalmass.get)
                if totalmass[best_clade] >= masscutoff:
                    monoclade[pref] = best_clade
                else:
                    monoclade[pref] = 'clade_not_conserved'
            else:
                print mode+'  UNRECOGNIZED MODE...(will be an error in future)'
    return monoclade

def write2fasta(ks_names, ks_seqs, outfile):
    f = open(outfile, "w")
    for i in range(0,len(ks_names),1):
        f.write(ks_names[i]+"\n"+ks_seqs[i]+"\n")
    f.close()
    return outfile
    
    

###LEGACY STUFF BELOW




# --------------------------------------
# -  core functions
# --------------------------------------
def alignment_per_ks(seq, name, data_dir, align_dir, nr):
    training_align = os.path.join(data_dir, "transKS_domains_annotated_alignedMuscle.fas")
    in_temp = os.path.join(align_dir, "in_seq" + str(nr) + ".fas")
    out_temp = os.path.join(align_dir, "out_seq" + str(nr) + ".fas")
    aligned = os.path.join(align_dir, "aligned" + str(nr) + ".fas")
    with open(in_temp, "w") as tmp_input:
        tmp_input.write("%s\n%s\n" % (name, seq))

    cmd = MuscleCommandline(profile='True', in1=training_align, in2=in_temp, out=out_temp)
    proc = subprocess.Popen(str(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    proc.communicate()
    retcode = proc.returncode

    if retcode == 1:
        print "Alignment failed"
    else:
        f_temp_input = open(out_temp, 'r').read()
        reformat(input=f_temp_input, out_filename=aligned)

    for f in [in_temp, out_temp]:
        os.remove(f)

    return aligned


# ----------------------------------------
# -- function: Phylotree_construct_perKS:
# ----------------------------------------
def phylotree_construct_perks_parallel(align_dir, tree_dir, max_core):
    fasfile_list = [os.path.join(align_dir, f) for f in os.listdir(align_dir) if
                    os.path.isfile(os.path.join(align_dir, f)) and not f.startswith('.')]

    out_tree_names = []
    cmds = []
    for f in fasfile_list:
        fname = os.path.split(f)[1]
        out_tree = os.path.join(tree_dir, ''.join((fname.split('.')[0], '.nwk')))
        cmd = 'fasttree %s > %s' % (f, out_tree)
        cmds.append(cmd)

    processes = (subprocess.Popen(c, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) for c in cmds)
    running_processes = list(islice(processes, max_core))  # start new processes

    logging.info("TransATPKS: constructing phylogeny tree for each KS domain")
    while running_processes:
        for i, process in enumerate(running_processes):
            err = process.communicate()[1]
            retcode = process.returncode
            f_name = err.split('\n')[1].split("/")[-1]
            tree_name = os.path.join(tree_dir, ''.join((f_name.split('.')[0], '.nwk')))
            if retcode == 0:
                out_tree_names.append(tree_name)
            else:
                logging.debug('%s:\n%s' % (f_name.split('.')[0], err))

            if process.poll() is not None:  # the process has finished
                try:
                    running_processes[i] = next(processes)  # start new process
                except StopIteration:
                    del running_processes[i]
                    break
    return out_tree_names


# ----------------------------------------
# -- function: add_feature_PhyloTree
# ----------------------------------------
def add_feature_PhyloTree(tree_input, data_dir):
    CladeInfo = os.path.join(data_dir, "all_KSinCompound_CladeInfo_v3_exclude_fasttree.txt")

    funClades = {}
    clade_info = open(CladeInfo).read().split('\n')[1:]
    for domain in clade_info:
        domain_id = domain.split('\t')[1]
        domain_id1 = re.sub(r'(\?)', '', domain_id)
        clade_id = domain.split('\t')[5]
        # clean up clades like II_c5_g1, but not out_group_like
        # So clean up things that look like roman numerals
        first_part = clade_id.split('_', 1)[0]
        if len(set(first_part) - set("IVXL")) == 0:
            clade_id = '_'.join(clade_id.split('_')[:2])
        funClades[domain_id1] = clade_id

    clade_list = [c for c in list(set(funClades.values())) if c != 'noClade']
    clade_length = []
    for c in clade_list:
        leaf_count = len([leaf_id for leaf_id, clade_id in funClades.items() if clade_id == c])
        clade_length.append(leaf_count)

    # Tree handel: add the new KS_id to dictionary funClades
    t = Tree(tree_input)
    leaf_id = [leaf.name for leaf in t]
    new_leaf = [l for l in leaf_id if l not in funClades.keys()]
    new_leaf_id = []

    for l in new_leaf:
        funClades[l] = 'new_%s' % (l.split('|')[-1])
        new_leaf_id.append('new_%s' % (l.split('|')[-1]))

    # add features, which are stored as value in dictionary funClades, to leaves of the tree
    for leaf in t:
        leaf.add_features(funClade=funClades.get(leaf.name, "none"))

    # add index to internal nodes
    edge = 1
    for node in t.traverse():
        if not node.is_leaf():
            node.name = 'Node_%d' % edge
            edge += 1

    return t, new_leaf_id, zip(clade_list, clade_length)


# --------------------------------------------
# -- Function: Score_Clade_Consistency
# --------------------------------------------
def Score_Clade_Consistency(tree, clade_summary, cutoff1, cutoff2):
    """
    #sample_input:
    tree_input = '/Users/Xiaowen/Desktop/tree/ABIB01000015_c1_KS_way1_part.nwk'
    clade_info = '/Users/Xiaowen/Documents/trans_AT/trans_AT_LuX_version2/all_KSinCompound_CladeInfo_v3_exclude.txt'
    t_feature, new_leaf , clade_summary = add_feature_PhyloTree(tree_input, CladeInfo)
    tree = t_feature
    clade_summary = clade_summary
    cutoff1 = 0.90
    cutoff2 = 0.90

    """
    clade_score = {}

    for clade in clade_summary:
        # get the nodes, whose leaves have feature as the clade id
        clade_id = clade[0]
        clade_len = clade[1]
        c1 = round(cutoff1 * clade_len)

        # get the node which only contain leaf belonging to the clade_id
        clade_i = tree.get_monophyletic(values=[clade_id], target_attr='funClade')

        clade_i_node = {}
        for node in clade_i:
            clade_i_node[node.name] = len(node.get_leaf_names())

        clade_size = clade_i_node.values()
        node_id_max = [id for id, node_len in clade_i_node.items() if node_len == max(clade_size)][0]

        if len(clade_size) == 1:
            recall = 1
            precision = 1
            clade_score[clade_id] = [recall, precision, node_id_max]
        elif len(clade_size) == 2:
            if max(clade_size) >= c1:
                recall = float(max(clade_size)) / clade_len
                precision = 1
                clade_score[clade_id] = [recall, precision, node_id_max]
            else:
                ancestor = tree.get_common_ancestor(clade_i_node.keys())
                recall = 1
                precision = float(clade_len) / len(ancestor.get_leaf_names())
                if precision >= cutoff2:
                    clade_score[clade_id] = [recall, precision, ancestor.name]
                else:
                    recall = float(max(clade_size)) / clade_len
                    clade_score[clade_id] = [recall, 1, node_id_max]
        elif len(clade_size) >= 3:
            if max(clade_size) >= c1:
                recall = float(max(clade_size)) / clade_len
                precision = 1
                clade_score[clade_id] = [recall, precision, node_id_max]
            else:
                node_id = [id for id in clade_i_node.keys() if id != node_id_max]
                dist = [tree.get_distance(node_id_max, i) for i in node_id]
                node_id_sorted = sort_XonY(node_id, dist, reverse=False)
                i = 0
                precision = 1
                while precision >= cutoff2 and i <= len(dist) + 1:
                    i += 1
                    ancestor = tree.get_common_ancestor([node_id_max] + node_id_sorted[0:i])
                    leaf_nr = sum([clade_i_node[n] for n in [node_id_max] + node_id_sorted[0:i]])
                    precision = float(leaf_nr) / len(ancestor.get_leaf_names())

                i_output = i - 1
                if i_output == 0:
                    recall = float(max(clade_size)) / clade_len
                    precision = 1
                    node_id_out = node_id_max
                else:
                    ancestor = tree.get_common_ancestor([node_id_max] + node_id_sorted[0:i_output])
                    leaf_nr = sum([clade_i_node[n] for n in [node_id_max] + node_id_sorted[0:i_output]])
                    precision = float(leaf_nr) / len(ancestor.get_leaf_names())
                    recall = float(leaf_nr) / clade_len
                    node_id_out = ancestor.name
                clade_score[clade_id] = [recall, precision, node_id_out]

    return clade_score


# ----------------------------------------------------
# -- function: assign_Functional_Clade_per_newDomain
# ----------------------------------------------------
def assign_Functional_Clade_per_newDomain(tree, new_leaf_ID, clade_summary,
                                          cutoff1, cutoff2,
                                          tree_conserve1, tree_conserve2,
                                          cutoff1_new, cutoff2_new):
    """
    #sample input:
    tree_input = ''
    tree_input = '/Users/Xiaowen/Documents/trans_AT/KS_domains_PredtransAT_treeML/ACOJ01000001_c13_KS_mixed_way3_part1.nwk'
    clade_info = '/Users/Xiaowen/Documents/trans_AT/trans_AT_LuX_version2/all_KSinCompound_CladeInfo_v3_exclude_fasttree.txt'
    t_feature, new_leaf , clade_summary = add_feature_PhyloTree(tree_input = tree_input, CladeInfo = clade_info)
    tree = t_feature
    print tree.get_ascii(attributes = ["name", "funClade"], show_internal = False) ##check the feature are add correctly
    cutoff1 = 0.90
    cutoff2 = 0.90
    new_leaf_ID = new_leaf[0]
    clade_info = clade_info
    tree_conserve1 = 0.8
    tree_conserve2 = 0.8
    cutoff2_new = 0.6
    cutoff1_new = 0.85

    """

    ori_clade_score = Score_Clade_Consistency(tree=tree, clade_summary=clade_summary, cutoff1=cutoff1, cutoff2=cutoff2)
    cutoff_r = [id for id, value in ori_clade_score.items() if value[0] >= cutoff1]
    cutoff_p = [id for id, value in ori_clade_score.items() if value[1] >= cutoff2]
    cutoff_r_p = list(set(cutoff_r) & set(cutoff_p))
    all_leaf = sum([i[1] for i in clade_summary])
    cutoff_r_p_leaf = sum(i[1] for i in clade_summary if i[0] in cutoff_r_p)
    conserve1 = float(len(cutoff_r_p)) / len(
        clade_summary)  # check whether most functional clades defined in original phylotree are conserved
    conserve2 = float(cutoff_r_p_leaf) / all_leaf  # check whether most leaves in original phylotree are conserved

    
    node_newKS = [node.name for node in tree.get_monophyletic(values=[new_leaf_ID], target_attr='funClade')][0]
    #print node_newKS + ' ' + node.name

    if conserve1 >= tree_conserve1 and conserve2 >= tree_conserve2:
        funClade_candidate = {}
        for clade in clade_summary:
            # print clade
            clade_id = clade[0]
            clade_len = clade[1]

            clade_i_score = ori_clade_score[clade_id]
            clade_i_leaf = [node.get_leaf_names() for node in tree.search_nodes(name=clade_i_score[2])][0]

            if node_newKS in clade_i_leaf:
                funClade_candidate[clade_id] = clade_i_score + [clade_len]
                break
            else:
                ancestor = tree.get_common_ancestor(clade_i_score[2], node_newKS)
                ancestor_funClade = [node.funClade for node in ancestor.get_leaves()]
                count1 = [c for c in ancestor_funClade if c == clade_id]
                recall = float(len(count1)) / clade_len
                precision = float(len(count1) + 1) / len(ancestor_funClade)
                funClade_candidate[clade_id] = [recall, precision, ancestor.name, len(ancestor_funClade)]

        funClade_candidate_filter = {}
        sum_score = []
        for k in funClade_candidate.keys():
            value_new = funClade_candidate[k]
            recall_ori = ori_clade_score[k][0]
            precision_ori = ori_clade_score[k][1]
            if value_new[0] >= min([cutoff1_new, recall_ori]) and value_new[1] >= min([cutoff2_new, precision_ori]):
                funClade_candidate_filter[k] = value_new + [value_new[0] + value_new[1]]
                sum_score.append(value_new[0] + value_new[1])

        if len(funClade_candidate_filter.keys()) >= 1:
            funClade_newKS_id = [id for id, value in funClade_candidate_filter.items() if value[-1] == max(sum_score)]
            funClade_newKS = '|'.join(funClade_newKS_id)
            recall_new = '|'.join([str(funClade_candidate[i][0]) for i in funClade_newKS_id])
            precision_new = '|'.join([str(funClade_candidate[i][1]) for i in funClade_newKS_id])
            clade_size_new = '|'.join([str(funClade_candidate[i][3]) for i in funClade_newKS_id])
        else:
            funClade_newKS = 'clade_not_conserved'
            recall_new = 'clade_not_conserved'
            precision_new = 'clade_not_conserved'
            clade_size_new = 'clade_not_conserved'
    else:
        funClade_newKS = 'tree_not_conserved'
        recall_new = 'tree_not_conserved'
        precision_new = 'tree_not_conserved'
        clade_size_new = 'tree_not_conserved'

    return node_newKS, funClade_newKS, conserve1, conserve2, recall_new, precision_new, clade_size_new


# -----------------------------------------
# --function: run_pipeline_perKS
# -----------------------------------------
def run_pipeline_perks_parallel(bgc_id, out_dir, align_dir, tree_dir, data_dir, cutoff1, cutoff2, tree_conserve1,
                                tree_conserve2, cutoff1_new, cutoff2_new, max_core):
    out_filename = os.path.join(out_dir, ''.join((bgc_id, '_funcAnnotate_perKS.txt')))
    out = open(out_filename, 'w')
    trees = phylotree_construct_perks_parallel(align_dir=align_dir, tree_dir=tree_dir, max_core=max_core)

    KS_func = {}
    if len(trees) >= 1:
        logging.info("Phylogenetic analysis: predicting substrate specificity of KS")
        for t in trees:
            t_feature, new_leaf, clade_summary = add_feature_PhyloTree(tree_input=t, data_dir=data_dir)
            KS_id, KS_clade, topo_conserve1, topo_conserve2, recall_new, precision_new, clade_size_new = assign_Functional_Clade_per_newDomain(
                tree=t_feature, new_leaf_ID=new_leaf[0], clade_summary=clade_summary, cutoff1=cutoff1, cutoff2=cutoff2,
                tree_conserve1=tree_conserve1, tree_conserve2=tree_conserve2, cutoff1_new=cutoff1_new,
                cutoff2_new=cutoff2_new)
            KS_func[KS_id] = [KS_clade, topo_conserve1, topo_conserve2, recall_new, precision_new, clade_size_new]

        for k in KS_func.keys():
            out.write('%s: %s\n' % (k, KS_func[k][0]))

        out.write('\n\n%s\n%s\t%s\t%s\t%s\t%s\t%s\n' % (
        '=' * 100, 'KS_id', 'Pencentage_conserved_clades', 'coverage_of_leaves', 'clade_recall', 'clade_precision',
        'clade_size'))

        for k in KS_func.keys():
            out.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (
            k, KS_func[k][1], KS_func[k][2], KS_func[k][3], KS_func[k][4], KS_func[k][5]))

        logging.debug('%s: successful for Phylogenetic analysis', bgc_id)

    else:
        logging.debug('%s: no tree is constructed', bgc_id)

    out.close()

    return KS_func

def run_pipeline_pplacer(bgc_id, out_dir, reference_alignment, alignment_file, data_dir, masscutoff, rt):
    out_filename = os.path.join(out_dir, ''.join((bgc_id, '_funcAnnotate_perKS.txt')))
    out = open(out_filename, 'w')
    
    logging.info("Phylogenetic analysis: predicting substrate specificity of KS")
    pplacer_tree = run_pplacer(reference_alignment, alignment_file, data_dir, rt)
    callmode = 'additive'
    clade_assignment = parse_pplacer(pplacer_tree, data_dir, masscutoff, callmode)
    for q in clade_assignment:
        out.write('%s: %s\n' % (q, clade_assignment[q]))
    logging.debug('%s: successful for Phylogenetic analysis', bgc_id)
    out.close()
    return clade_assignment

ks_names = []
ks_seqs = []
fa = SeqIO.parse(open(sys.argv[1]),'fasta')
rt = sys.argv[2]
d = './data'
ra = sys.argv[3]
for rec in fa:
    ks_names.append('>'+rec.id) ## NOTE: this assumes this header will be unique when compared to the reference set. if not, will return OUTGROUP regardless of the call. If running such non-uniue seqs, change the header to '>query_'+rec.id or similar make it unique
    ks_seqs.append(rec.seq)
    alignment_file = align_ks_domains(ra, ks_names, ks_seqs, d)
    clade = run_pipeline_pplacer('precomp', './', ra, alignment_file, d, 0.6, rt)
    for k in clade:
        print "\t".join([str(k), str(clade[k])])
