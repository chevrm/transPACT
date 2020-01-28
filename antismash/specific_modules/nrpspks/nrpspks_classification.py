import os
import re
from os import path
from helperlibs.wrappers.io import TemporaryDirectory
from antismash import utils
from orderfinder import find_clusterpksnrpsgenes
from nrpspksdomainalign import substrates_specificity_phylo as ssp
from nrpspksdomainalign import index_domains as index_domains
from nrpspksdomainalign import get_msa_dist_matrix as mdm
from nrpspksdomainalign import calculate_distance_transATPKS as cdt
from nrpspksdomainalign import plot_domain_align as pda
from datetime import datetime

def _get_transatpks_geneclusters(pksnrpsvars, seq_record):
    nrpspksclusters = list(set(utils.get_cluster_features_of_type(seq_record, "transatpks")))
    genes_in_cluster = {}
    for cluster in nrpspksclusters:
        cluster_id = utils.get_cluster_number(cluster)
        cluster_genes = [utils.get_gene_id(feature) for feature in find_clusterpksnrpsgenes(cluster, pksnrpsvars.pksnrpscoregenes)]
        genes_in_cluster[cluster_id] = cluster_genes
    return genes_in_cluster

def _get_nrpspks_domains_ks(pksnrpsvars, seq_record, domain):
    transatpks_geneclusters = _get_transatpks_geneclusters(pksnrpsvars, seq_record)
    transatpks_genes = list(set([g for g_list in transatpks_geneclusters.values() for g in g_list]))
    ksnames = []
    ksseqs = []
    if len(transatpks_geneclusters) >= 1:
        job_id = seq_record.id
        for feature in pksnrpsvars.pksnrpscoregenes:
            start_cds = str(feature.location.nofuzzy_start)
            end_cds = str(feature.location.nofuzzy_end)
            strand = feature.location.strand
            if strand == 1:
                strand_char = '+'
            else:
                strand_char = '-'
            loc = '-'.join((start_cds, end_cds))
            prot_id = product = ''
            if 'protein_id' in feature.qualifiers:
                prot_id = feature.qualifiers["protein_id"][0]
            if 'product' in feature.qualifiers:
                product = feature.qualifiers["product"][0].replace(' ', '_').replace('|', '')
                # We use | as a separator later
                assert '|' not in product, product
            gene_id = utils.get_gene_id(feature)
            if gene_id in transatpks_genes:
                domaindetails = pksnrpsvars.domaindict[gene_id]
                nr = 0
                for tab in domaindetails:
                    if tab[0] == domain:
                        nr += 1
                        start = int(tab[1])
                        end = int(tab[2])
                        loc_domain = '-'.join((str(start), str(end)))
                        ks_index = ''.join(('KS', str(nr)))
                        name1 = '|'.join(
                            [''.join(['>', job_id]), 'c', loc, strand_char, gene_id, product, prot_id, loc_domain, ks_index])
                        name = re.sub(r'(\:|\'|\(|\)|\,|\?|\;)', '', name1)
                        seq = str(utils.get_aa_sequence(feature))[start:end]
                        ksnames.append(name)
                        ksseqs.append(seq)
    return ksnames, ksseqs

def _index_domains_per_cluster(genes_per_cluster, domains_annotation, domain_ab):
    domains_per_cluster = {}
    for c in genes_per_cluster.keys():
        k = c
        v = {}
        gene_list = genes_per_cluster[c]
        for d in domains_annotation.keys():
            d_info = d.split("|")
            d_info[1] = d_info[1]+str(k)
            new_d = '|'.join(d_info)
            g = d.split('|')[4]
            if g in gene_list:
                v[new_d] = domains_annotation[d]
        v_domain_index = index_domains.run_index_domain(domain_id_list=v.keys(), domain_ab = domain_ab)
        v_new = {}
        for id_old in v.keys():
            id_new = v_domain_index[id_old]
            v_new[id_new] = v[id_old]
        domains_per_cluster[k] = v_new
    return domains_per_cluster

def classify_nrpspks_domains_ks(pksnrpsvars, seq_record, options):
    print str(datetime.now())+"\tBegin classify_nrpspks_domains_ks"
    #todo: add an argument for the number of selected BGCs in run_antismash via parser.add_argument(),the parsed info will be store in options
    bgcs_nr = options.transatpks_da_cutoff
    with TemporaryDirectory(change=True):
        options.classify_domain_outputfolder_align = path.abspath(path.join(options.raw_predictions_outputfolder, "classified_domain", "align_ks"))
        if not os.path.exists(options.classify_domain_outputfolder_align):
            os.makedirs(options.classify_domain_outputfolder_align)
        options.classify_domain_outputfolder_tree = path.abspath(path.join(options.raw_predictions_outputfolder, "classified_domain", "tree_ks"))
        if not os.path.exists(options.classify_domain_outputfolder_tree):
            os.makedirs(options.classify_domain_outputfolder_tree)
        ks_names, ks_seqs = _get_nrpspks_domains_ks(pksnrpsvars, seq_record, domain="PKS_KS")   #todo: write ks_names and ks_seqs into a fasta file
        genes_per_transatpks = _get_transatpks_geneclusters(pksnrpsvars, seq_record)
        nrpspksdomainalign_dir = utils.get_full_path(__file__, "nrpspksdomainalign")
        data_dir = path.join(nrpspksdomainalign_dir, "data")
        annotation_out_dir = os.path.split(options.classify_domain_outputfolder_align)[0]
        genomeSize_info_seq_record = len(seq_record.seq)
        speciesName_seq_record = 'Query sequence'
        #Step 1: align new sequences to reference alignment, also return distance matrix
        print str(datetime.now())+"\tAlign to reference"
        reference_alignment = os.path.join(data_dir, "647KS_mcformat.afa")
        alignment_file = ssp.align_ks_domains(reference_alignment, ks_names, ks_seqs, data_dir)
        #Step 2: Run pplacer
        #Step 3: Assign clades using ETE
        print str(datetime.now())+"\tAssign clades"
        KS_annotation = ssp.run_pipeline_pplacer(bgc_id = seq_record.id, 
                                                 out_dir = annotation_out_dir, 
                                                 reference_alignment = reference_alignment,
                                                 alignment_file = alignment_file,
                                                 data_dir = data_dir,
                                                 masscutoff = 0.9)
        #for k in KS_annotation:
        #    print "\t".join([str(k), str(KS_annotation[k])])
        #Step 4: Find nearest neighbours in reference assembly lines
        print str(datetime.now())+"\tCompute nearest neighbor matrix"
        KS_per_cluster = _index_domains_per_cluster(genes_per_transatpks, KS_annotation, "KSA")
        precomputed_aln = os.path.join(data_dir, "KS_precomputed_hmmalign.afa")
        precomputed_dist = os.path.join(data_dir, "precomputed_fracid.csv")
        KS_msa_dist = mdm.run_get_msa_dist_matrix_precomputed(ks_names, ks_seqs, data_dir,"KS_rawseq_pred_training_transATPKS.txt", annotation_out_dir, KS_per_cluster, precomputed_aln, precomputed_dist)
        #Step 5: Align query with nearest neighbours (local distance all vs all)
        print str(datetime.now())+"\tCalculate distance"
        similar_bgc_per_cluster, new_cluster, new_cluster_index = cdt.run_calculate_distance(data_dir, KS_msa_dist, KS_per_cluster, bgcs_nr)
        #Step 6: Visualize results
        print str(datetime.now())+"\tVisualize results"
        for b in similar_bgc_per_cluster.keys():
            b_id = re.sub("c",'', b.split('|')[1])
            dist_df = similar_bgc_per_cluster[b][1]
            #annotateKS_training = os.path.join(data_dir, "annotateKS_perCompound_v4_uniPK.txt")
            annotateKS_training = os.path.join(data_dir, "annotateKS_training.txt")
            #annotateKS_novel = os.path.join(data_dir, "annotateKS_per_pred_transATPKS.txt")
            annotateKS_novel = os.path.join(data_dir, "annotateKS_novel.txt")
            #color_code_file = os.path.join(data_dir, "Clade_PhyloNode_v3_color.txt")
            color_code_file = os.path.join(data_dir, "transAT_mc.txt")
            #KS_info_file_list = [os.path.join(data_dir, file) for file in ['KSinCompound_CladeInfo_v4_exclude.txt', 'KSinCompound_CladeInfo_pred_transATPKS.txt']]
            KS_info_file_list = [os.path.join(data_dir, file) for file in ['KS_fullclades.infolist.tsv']]
            genomeSize_file_list = [os.path.join(data_dir, file) for file in['genome_size_training_transATPKS.txt', 'genome_size_pred_transATPKS.txt']]
            speciesName_file_list = [os.path.join(data_dir, file) for file in['species_ID_training_transATPKS.txt', 'species_ID_pred_transATPKS.txt']]
            outfile_svg_temp = os.path.join(annotation_out_dir, ''.join(["domain_align_cluster", b_id, ".svg"]))
            svg_fout = os.path.join(annotation_out_dir, ''.join(["domain_align_cluster", b_id, "_modify.svg"]))
            pda.run_plot_network_BGC_list(dist_df=dist_df,
                                  new_bgc_id=[b],
                                  new_cluster=new_cluster,
                                  new_cluster_index=new_cluster_index,
                                  annotateKS_training=annotateKS_training,
                                  annotateKS_novel=annotateKS_novel,
                                  speciesName_file_list=speciesName_file_list,
                                  dist=0.2,
                                  color_code_file=color_code_file,
                                  outfile=outfile_svg_temp,
                                  plt_width=40,
                                  KS_info_file_list=KS_info_file_list,
                                  genomeSize_file_list=genomeSize_file_list,
                                  ksa_per_new_cluster=KS_per_cluster,
                                  genomeSize_new_bgc=genomeSize_info_seq_record ,
                                  speciesName_new_bgc=speciesName_seq_record,
                                  seq_simScore = KS_msa_dist,
                                  svg_fout=svg_fout)
        print str(datetime.now())+"\tDONE!"
        import sys
        sys.exit()
        return
