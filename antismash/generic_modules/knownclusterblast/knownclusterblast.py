# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2010-2014 Marnix H. Medema
# University of Groningen
# Department of Microbial Physiology / Groningen Bioinformatics Centre
#
# Copyright (C) 2011,2012 Kai Blin
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Div. of Microbiology/Biotechnology
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
import os
from antismash import utils
from ..clusterblast import clusterblast
from helperlibs.wrappers.io import TemporaryDirectory

def perform_knownclusterblast(options, seq_record, clusters, proteins):
    # Run BLAST on gene cluster proteins of each cluster and parse output
    logging.debug("Running DIAMOND knowncluster searches..")
    geneclusters = utils.get_sorted_cluster_features(seq_record)

    all_names, all_seqs, all_prots = [], [] ,[]
    prots_by_cluster = []
    for genecluster in geneclusters:
        names, seqs, prots = clusterblast.create_blast_inputs(genecluster, seq_record)
        all_names.extend(names)
        all_seqs.extend(seqs)
        all_prots.extend(prots)
        prots_by_cluster.append(prots)

    debug_path = os.path.join(options.dbgclusterblast, "knownclusterblastoutput.txt")
    if options.dbgclusterblast and os.path.exists(debug_path):
        logging.debug("Skipping DIAMOND calculations, using previous results")
        with open(debug_path, "r") as fh:
            blastoutput = fh.read()
    else:
        with TemporaryDirectory(change=True) as tempdir:
            utils.writefasta([qcname.replace(" ", "_") for qcname in all_names],
                             all_seqs, "input.fasta")
            out, err, retcode = clusterblast.run_diamond("input.fasta",
                                            os.path.join(options.knownclusterblastdir, 'knownclusterprots'),
                                            tempdir, options)
            if retcode != 0:
                logging.debug("out: %r, err: %r, retcode: %s", out, err, retcode)
            with open("input.out", 'r') as fh:
                blastoutput = fh.read()
            clusterblast.write_raw_clusterblastoutput(options.full_outputfolder_path, blastoutput,
                                         searchtype="knownclusters")

    minseqcoverage = 40
    minpercidentity = 45
    clusters_by_number, _ = clusterblast.parse_all_clusters(blastoutput, minseqcoverage,
                                              minpercidentity, seq_record)

    knownclusterblastStorage = utils.Storage()
    knownclusterblastStorage.clusters = clusters
    knownclusterblastStorage.proteins = proteins

    for genecluster, queryclusterprots in zip(geneclusters, prots_by_cluster):
        clusternumber = utils.get_cluster_number(genecluster)
        cluster_names_to_queries = clusters_by_number.get(clusternumber, {})
        allcoregenes = [utils.get_gene_id(cds) for cds in utils.get_secmet_cds_features(seq_record)]
        ranking = clusterblast.score_clusterblast_output(clusters, allcoregenes, cluster_names_to_queries)

        # store all clusterblast related data in a utils.Storage object and serialize it
        knownclusterblastStorage.clusternumber = clusternumber
        knownclusterblastStorage.queryclusterprots = queryclusterprots
        knownclusterblastStorage.ranking = ranking
        clusterblast.write_clusterblast_output(options, seq_record, knownclusterblastStorage, searchtype="knownclusters")

    mibig_protein_homology(blastoutput, seq_record, geneclusters, clusters, options)


def mibig_protein_homology(blastoutput, seq_record, geneclusters, clusters, options):

    minseqcoverage = 20
    minpercidentity = 20
    _, queries_by_cluster = clusterblast.parse_all_clusters(blastoutput, minseqcoverage,
                                               minpercidentity, seq_record)

    for genecluster in geneclusters:
        cluster_number = utils.get_cluster_number(genecluster)
        queries = queries_by_cluster.get(cluster_number, {})

        # Since the BLAST query was only for proteins in the cluster just need to iterate through the keys and generate
        # a file for each of the keys
        outputfolder = os.path.join(options.knownclusterblast_outputfolder,
                                    "cluster{}".format(cluster_number))
        if not os.path.exists(outputfolder):
            os.mkdir(outputfolder)
        for cluster_protein in queries.values():
            protein_name = cluster_protein.id
            with open(outputfolder + os.sep + protein_name + '_mibig_hits.txt', 'w') as outfile:
                outfile.write('#Protein\tDescription\tMiBIG Cluster\tMiBIG Product'
                              '\tPercent ID\tPercent Coverage\tBLAST Score\t Evalue\n')
                for subject in cluster_protein.subjects.values():
                    gene_id = subject.locus_tag
                    gene_descr = subject.annotation
                    mibig_cluster = subject.genecluster
                    mibig_product = clusters[mibig_cluster][1]
                    percent_id = str(subject.perc_ident)
                    blast_score = str(subject.blastscore)
                    percent_cvg = str(subject.perc_coverage)
                    e_value = str(subject.evalue)
                    outfile.write(gene_id + '\t' + gene_descr + '\t' + mibig_cluster
                                  + '\t' + mibig_product + '\t' + percent_id + '\t' + percent_cvg
                                  + '\t' + blast_score + '\t' + e_value + '\n')
