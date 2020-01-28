# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2010-2012 Marnix H. Medema
# University of Groningen
# Department of Microbial Physiology / Groningen Bioinformatics Centre
#clusterblastStorage = utils.Storage()
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

def perform_subclusterblast(options, seq_record, clusters, proteins):
    #Run BLAST on gene cluster proteins of each cluster and parse output
    logging.debug("Running NCBI BLAST+ subcluster searches..")
    geneclusters = utils.get_sorted_cluster_features(seq_record)
    with TemporaryDirectory(change=True):
        for genecluster in geneclusters:
            clusternumber = utils.get_cluster_number(genecluster)
            debug_path = os.path.join(options.dbgclusterblast, "subclusterblast",
                                      "cluster" + str(clusternumber) + ".txt")
            if options.debug and os.path.exists(debug_path):
                logging.debug("Skipping SubClusterblast calculations, using results from %s instead", debug_path)
            else:
                logging.debug("   Gene cluster " + str(clusternumber))
                queryclusternames, queryclusterseqs, queryclusterprots = clusterblast.create_blast_inputs(genecluster, seq_record)
                clusterblast.write_clusterblast_inputfiles(options, queryclusternames, queryclusterseqs)
                clusterblast.run_clusterblast_processes(options, searchtype="subclusters")
                blastoutput = clusterblast.read_clusterblast_output(options)
                clusterblast.write_raw_clusterblastoutput(options.full_outputfolder_path, blastoutput, searchtype="subclusters")
                logging.debug("   Blast search finished. Parsing results...")
                minseqcoverage = 40
                minpercidentity = 45
                _, cluster_names_to_queries = clusterblast.blastparse(blastoutput, minseqcoverage, minpercidentity, seq_record)
                allcoregenes = [utils.get_gene_acc(cds) for cds in utils.get_secmet_cds_features(seq_record)]
                ranking = clusterblast.score_clusterblast_output(clusters, allcoregenes, cluster_names_to_queries)
                
                # store all clusterblast related data in a utils.Storage object and serialize it
                subclusterblastStorage = utils.Storage()
                subclusterblastStorage.clusternumber = clusternumber
                subclusterblastStorage.queryclusterprots = queryclusterprots
                subclusterblastStorage.clusters = clusters
                subclusterblastStorage.ranking = ranking
                subclusterblastStorage.proteins = proteins
                    
                clusterblast.write_clusterblast_output(options, seq_record, subclusterblastStorage, searchtype="subclusters")
