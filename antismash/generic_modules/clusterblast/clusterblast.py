# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2010-2012 Marnix H. Medema
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
import sys
import os
from os import path
from collections import OrderedDict
from antismash import utils
from multiprocessing import Process
from antismash import config
from helperlibs.wrappers.io import TemporaryDirectory

def runblast(query, target):
    command = ["blastp", "-db", target, "-query", query, "-outfmt", "6", "-max_target_seqs", "10000", "-evalue", "1e-05", "-out", query.split(".")[0] + ".out"]
    utils.execute(command)


def run_diamond(query, target, tempdir, options):
    command = [
        "diamond", "blastp",
        "--db", target,
        "--threads", str(options.cpus),
        "--query", query,
        "--compress", "0",
        "--max-target-seqs", "10000",
        "--evalue", "1e-05",
        "--out", "input.out",
        "--outfmt", "6",  # 6 is blast tabular format, just as in blastp
        "--tmpdir", tempdir
    ]
    return utils.execute(command)


def make_blastdb(inputfile, dbname):
    command = ["makeblastdb", "-in", inputfile, "-out", dbname, "-dbtype", "prot"]
    utils.execute(command)

def load_geneclusters(searchtype):
    #Load gene cluster database into memory
    options = config.get_config()
    if not 'clusterblastdir' in options:
        options.clusterblastdir = path.dirname(utils.get_full_path(__file__, ''))
        options.subclusterblastdir = path.join(path.dirname(options.clusterblastdir), 'subclusterblast')
        options.knownclusterblastdir = path.join(path.dirname(options.clusterblastdir), 'knownclusterblast')
    else:
        options.subclusterblastdir = path.join(path.dirname(path.dirname(utils.get_full_path(__file__, ''))), 'subclusterblast')
        options.knownclusterblastdir = path.join(path.dirname(path.dirname(utils.get_full_path(__file__, ''))), 'knownclusterblast')

    if searchtype == "general":
        logging.debug("ClusterBlast: Loading gene clusters database into memory...")
        geneclustersfile = path.join(options.clusterblastdir, "geneclusters.txt")
    elif searchtype == "subclusters":
        logging.debug("SubClusterBlast: Loading gene clusters database into memory...")
        geneclustersfile = path.join(options.subclusterblastdir, "subclusters.txt")
    elif searchtype == "knownclusters":
        logging.debug("KnownClusterBlast: Loading gene clusters database into memory...")
        geneclustersfile = path.join(options.knownclusterblastdir, "knownclusters.txt")
    geneclustersfile = open(geneclustersfile,"r")
    filetext = geneclustersfile.read()
    lines = [line for line in filetext.split("\n") if "\t" in line]
    clusters = {}
    for i in lines:
        tabs = i.split("\t")
        accession = tabs[0]
        clusterdescription = tabs[1]
        clusternr = tabs[2]
        clustertype = tabs[3]
        clustername = accession + "_" + clusternr
        clustertags = tabs[4].split(";")
        clusterprots = tabs[5].split(";")
        clusters[clustername] = [clusterprots,clusterdescription,clustertype,clustertags]
    return clusters

class Protein(object):
    """ Holds details of a protein. Members cannot be added dynamically.
    """
    # At time of writing this class will be instantiated ~7 million times per
    # clusterblast invocation.
    # With those numbers, the memory use of the class without slots: 2.3 Gb
    #                                                and with slots: 0.4 Gb
    __slots__ = ("name", "locus_tag", "location", "strand", "annotations")
    def __init__(self, name, locus_tag, location, strand, annotations):
        self.name = name
        self.locus_tag = locus_tag
        self.location = location
        self.strand = strand
        self.annotations = annotations

    def __str__(self):
        if len(self.location.split("-")) != 2:
            raise ValueError("Invalid location in Protein: %s"%self.location)
        tag = self.locus_tag
        if tag == "no_locus_tag":
            tag = self.name
        locations = self.location.replace("-", "\t")
        return "{}\t{}\t{}\t{}\t{}\n".format(tag, self.name, locations,
                                                 self.strand, self.annotations)

def load_geneclusterproteins(accessiondict, searchtype):
    options = config.get_config()
    if not 'clusterblastdir' in options:
        options.clusterblastdir = path.dirname(utils.get_full_path(__file__, ''))
        options.subclusterblastdir = path.join(path.dirname(options.clusterblastdir), 'subclusterblast')
        options.knownclusterblastdir = path.join(path.dirname(options.clusterblastdir), 'knownclusterblast')
    else:
        options.subclusterblastdir = path.join(path.dirname(path.dirname(utils.get_full_path(__file__, ''))), 'subclusterblast')
        options.knownclusterblastdir = path.join(path.dirname(path.dirname(utils.get_full_path(__file__, ''))), 'knownclusterblast')
    #Load gene cluster database proteins info into memory
    if searchtype == "general":
        logging.debug("ClusterBlast: Loading gene cluster database proteins into " \
        "memory...")
        gclusterprotsfile = path.join(options.clusterblastdir, "geneclusterprots.fasta")
    elif searchtype == "subclusters":
        logging.debug("SubClusterBlast: Loading gene cluster database proteins into " \
        "memory...")
        gclusterprotsfile = path.join(options.subclusterblastdir, "subclusterprots.fasta")
    elif searchtype == "knownclusters":
        logging.debug("KnownClusterBlast: Loading gene cluster database proteins into " \
        "memory...")
        gclusterprotsfile = path.join(options.knownclusterblastdir, "knownclusterprots.fasta")

    proteins = {}

    with open(gclusterprotsfile, 'r') as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line[0] != ">":
                continue
            tabs = line.split("|")
            locustag = tabs[4]
            if accessiondict.has_key(locustag):
                locustag = "h_" + locustag
            location = tabs[2]
            strand = tabs[3]
            annotations = tabs[5]
            name = tabs[6]
            proteins[name] = Protein(name, locustag, location, strand, annotations)
    return proteins

def load_clusterblast_database(seq_record, searchtype="general"):
    accessiondict = {}
    for cds in utils.get_cds_features(seq_record):
        accessiondict[utils.get_gene_acc(cds)] = utils.get_gene_accession(cds)
    clusters = load_geneclusters(searchtype)
    proteins = load_geneclusterproteins(accessiondict, searchtype)
    return clusters, proteins

def create_blast_inputs(genecluster, seq_record):
    #Create input fasta files for BLAST search
    queryclusterprots = utils.get_cluster_cds_features(genecluster, seq_record)
    queryclusternames = []
    queryclusterseqs = []
    queryclusterprotsnames = []
    for cds in queryclusterprots:
        if cds.strand == 1:
            strand = "+"
        else:
            strand = "-"
        fullname = "|".join(["input", "c" + str(utils.get_cluster_number(genecluster)), \
                             str(cds.location.nofuzzy_start) + "-" + \
                             str(cds.location.nofuzzy_end), \
                             strand, utils.get_gene_acc(cds), utils.get_gene_annotation(cds)])
        queryclusterseqs.append(str(utils.get_aa_sequence(cds)))
        queryclusternames.append(fullname)
        queryclusterprotsnames.append(utils.get_gene_acc(cds))

    return queryclusternames, queryclusterseqs, queryclusterprotsnames

def run_internal_blastsearch():
    #Run and parse BLAST search
    make_blastdb("internal_input.fasta", "internal_input.fasta")
    runblast("internal_input.fasta", "internal_input.fasta")
    blastoutput = open("internal_input.out","r").read()
    return blastoutput

def uniqueblasthitfilter(blastlines):
    #Filter for best blast hits (of one query on each subject)
    query_subject_combinations = set()
    blastlines2 = []
    for tabs in blastlines:
        # if it doesn't even have both values, it's not a hit so skip it
        if len(tabs) < 2:
            continue
        query = tabs[0]
        subject = tabs[1]
        query_subject_combination = (query, subject)
        if query_subject_combination not in query_subject_combinations:
            query_subject_combinations.add(query_subject_combination)
            blastlines2.append(tabs)
    return blastlines2

class Subject(object):
    def __init__(self, name, genecluster, start, end, strand, annotation,
                 perc_ident, blastscore, perc_coverage, evalue, locus_tag):
        self.name = name
        self.genecluster = genecluster
        self.start = start
        self.end = end
        self.strand = strand
        self.annotation = annotation
        self.perc_ident = int(perc_ident)
        self.blastscore = int(blastscore)
        self.perc_coverage = float(perc_coverage)
        self.evalue = float(evalue)
        self.locus_tag = locus_tag

    def get_table_string(self):
        return "\t".join([str(i) for i in [self.name, self.perc_ident,
                                           self.blastscore, self.perc_coverage,
                                           self.evalue]])

class Query(object):
    def __init__(self, entry, index):
        parts = entry.split("|")
        self.cluster_number = int(parts[1][1:]) # c1 -> 1
        self.id = parts[4]
        self.entry = entry
        self.subjects = OrderedDict()
        self.cluster_name_to_subjects = {}
        self.index = index

    def add_subject(self, subject):
        self.subjects[subject.name] = subject
        if subject.genecluster not in self.cluster_name_to_subjects:
            self.cluster_name_to_subjects[subject.genecluster] = []
        self.cluster_name_to_subjects[subject.genecluster].append(subject)

    def get_subjects_by_cluster(self, cluster_name):
        return self.cluster_name_to_subjects.get(cluster_name, [])

def parse_subject(tabs, seqlengths, geneclustergenes, seq_record):
    if len(tabs) < 12:
        logging.error("Malformed blast pairing: %s", "\t".join(tabs))
    query = tabs[0]
    subject_parts = tabs[1].split("|")
    subject = subject_parts[4]
    if subject == "no_locus_tag":
        subject = subject_parts[6]
    if subject in geneclustergenes:
        subject = "h_" + subject
    if len(subject_parts) > 6:
        locustag = subject_parts[6]
    else:
        locustag = ""
    genecluster = "{}_{}".format(subject_parts[0], subject_parts[1])
    start, end = subject_parts[2].split("-")[:2]
    strand = subject_parts[3]
    annotation = subject_parts[5]
    perc_ident = int(float(tabs[2]) + 0.5)
    evalue = str(tabs[10])
    blastscore = int(float(tabs[11])+0.5)
    query_key = query.split("|")[4]
    if seqlengths.has_key(query_key):
        perc_coverage = (float(tabs[3]) / seqlengths[query_key]) * 100
    else:
        feature_by_id = utils.get_feature_dict_protein_id(seq_record)
        seqlength = len(utils.get_aa_sequence(feature_by_id[query_key]))
        perc_coverage = (float(tabs[3]) / seqlength) * 100
    return Subject(subject, genecluster, start, end, strand, annotation,
                   perc_ident, blastscore, perc_coverage, evalue, locustag)

def parse_all_clusters(blasttext, minseqcoverage, minpercidentity, seq_record):
    """ Parses blast results, groups into results by cluster number

        blasttext: the output from diamond in blast format
        minseqcoverage: the exclusive lower bound of sequence coverage for a match
        minpercidentity: the exclusive lower bound of identity similarity for a match
        seq_record: used to get all gene ids in the cluster, and used as a
                backup to fetch sequence length if missing from seqlengths
    """
    seqlengths = fastaseqlengths(seq_record)
    geneclustergenes = [utils.get_gene_acc(cds) for cds in utils.get_withincluster_cds_features(seq_record)]
    queries = OrderedDict()
    clusters = OrderedDict()
    blastlines = uniqueblasthitfilter([line.split("\t") for line in blasttext.rstrip().split("\n")])
    current_query = None
    queries_by_cluster_number = {}
    clusters_by_query_cluster_number = {}

    for tabs in blastlines:
        query = tabs[0]
        subject = parse_subject(tabs, seqlengths, geneclustergenes, seq_record)

        # only process the pairing if limits met
        if subject.perc_ident <= minpercidentity \
                or subject.perc_coverage <= minseqcoverage:
            continue

        new_query = query not in queries

        if new_query:
            current_query = Query(query, len(queries))
            cluster_number = current_query.cluster_number
            # is it a new cluster number? if so, reset collections
            if cluster_number not in queries_by_cluster_number:
                queries = OrderedDict()
                clusters = OrderedDict()
                # reset query index, since we started a new collection
                current_query.index = 0
                # link them
                queries_by_cluster_number[cluster_number] = queries
                clusters_by_query_cluster_number[cluster_number] = clusters
            # finally, add the query to the current tracker
            queries[query] = current_query

        if subject.genecluster not in clusters:
            clusters[subject.genecluster] = []
        clusters[subject.genecluster].append(current_query)

        # link the subject to the query
        current_query.add_subject(subject)

    return clusters_by_query_cluster_number, queries_by_cluster_number

def blastparse(blasttext, minseqcoverage, minpercidentity, seq_record):
    """ blasttext: the output from diamond in blast format
        minseqcoverage: the exclusive lower bound of sequence coverage for a match
        minpercidentity: the exclusive lower bound of identity similarity for a match
        seq_record: used to get all gene ids in the cluster, and used as a
                backup to fetch sequence length if missing from seqlengths
    """
    seqlengths = fastaseqlengths(seq_record)
    geneclustergenes = [utils.get_gene_acc(cds) for cds in utils.get_withincluster_cds_features(seq_record)]
    queries = OrderedDict()
    clusters = OrderedDict()
    blastlines = uniqueblasthitfilter([line.split("\t") for line in blasttext.rstrip().split("\n")])
    current_query = None

    for tabs in blastlines:
        query = tabs[0]
        subject = parse_subject(tabs, seqlengths, geneclustergenes, seq_record)

        # only process the pairing if limits met
        if subject.perc_ident <= minpercidentity \
                or subject.perc_coverage <= minseqcoverage:
            continue

        new_query = query not in queries
        new_hit = subject.genecluster not in clusters

        if new_query:
            current_query = Query(query, len(queries))
            queries[query] = current_query

        if new_hit:
            clusters[subject.genecluster] = []
        clusters[subject.genecluster].append(current_query)

        # link the subject to the query
        current_query.add_subject(subject)

    return queries, clusters

def fastaseqlengths(seq_record):
    seqlengths = {}
    cdsfeatures = utils.get_cds_features(seq_record)
    for cds in cdsfeatures:
        seqlength = len(str(utils.get_aa_sequence(cds)))
        seqlengths[utils.get_gene_acc(cds)] = seqlength
    return seqlengths

def find_internal_orthologous_groups(queries, clusternames):
    #find and store internal homologs
    groups = []
    for name in clusternames:
        if name not in queries:
            groups.append([name.split("|")[4]])
            continue
        query = queries[name]
        group = []
        for hit in query.subjects:
            if hit.startswith("h_"):
                group.append(hit[2:])
            else:
                group.append(hit)
        if query.id not in group:
            group.append(query.id)
        x = 0
        for l in groups:
            for m in group:
                if m in l:
                    del groups[x]
                    for n in l:
                        if n not in group:
                            group.append(n)
                    break
            x += 1
        group.sort()
        groups.append(group)
    return groups

def internal_homology_blast(seq_record):
    #Run BLAST on gene cluster proteins of each cluster on itself to find internal homologs, store groups of homologs - including singles - in a dictionary as a list of lists accordingly
    with TemporaryDirectory(change=True):
        logging.debug("Finding internal homologs in each gene cluster..")
        internalhomologygroups = {}
        geneclusters = utils.get_sorted_cluster_features(seq_record)
        for genecluster in geneclusters:
            clusternumber = utils.get_cluster_number(genecluster)
            iqueryclusternames, iqueryclusterseqs, _ = create_blast_inputs(genecluster, seq_record)
            utils.writefasta(iqueryclusternames, iqueryclusterseqs, "internal_input.fasta")
            blastoutput = run_internal_blastsearch()
            queries, _ = blastparse(blastoutput, 25, 30, seq_record)
            groups = find_internal_orthologous_groups(queries, iqueryclusternames)
            internalhomologygroups[clusternumber] = groups
    return internalhomologygroups

def write_clusterblast_inputfiles(options, queryclusternames, queryclusterseqs):
    equalpartsizes = int(len(queryclusternames) / options.cpus)
    for i in range(options.cpus):
        if i == 0:
            setnames = queryclusternames[:equalpartsizes]
            setseqs = queryclusterseqs[:equalpartsizes]
        elif i == (options.cpus - 1):
            setnames = queryclusternames[(i*equalpartsizes):]
            setseqs = queryclusterseqs[(i*equalpartsizes):]
        else:
            setnames = queryclusternames[(i*equalpartsizes):((i+1)*equalpartsizes)]
            setseqs = queryclusterseqs[(i*equalpartsizes):((i+1)*equalpartsizes)]
        utils.writefasta(setnames, setseqs, "input" + str(i) + ".fasta")

def run_clusterblast_processes(options, searchtype="general"):

    processes = []
    for i in range(options.cpus):
        if searchtype == "general":
            blastdb = path.join(options.clusterblastdir, 'geneclusterprots.fasta')
        elif searchtype == "subclusters":
            blastdb = path.join(options.subclusterblastdir, 'subclusterprots.fasta')
        elif searchtype == "knownclusters":
            blastdb = path.join(options.knownclusterblastdir, 'knownclusterprots.fasta')
        if sys.platform == 'win32':
            blastpath, _, blastdb = blastdb.rpartition(os.sep)
            os.environ['BLASTDB'] = blastpath
        processes.append(Process(target=runblast, args=["input" + str(i) + ".fasta", blastdb]))
    for i in processes:
        i.start()
    for i in processes:
        i.join()

def read_clusterblast_output(options):
    blastoutput = []
    for i in range(options.cpus):
        fh = open("input" + str(i) + ".out","r")
        output = fh.read()
        fh.close()
        blastoutput.append(output)
    return "".join(blastoutput)

def write_raw_clusterblastoutput(outputfoldername, blastoutput, searchtype="general"):
    if searchtype == "general":
        blastoutputfile = open(outputfoldername + os.sep + "clusterblastoutput.txt","w")
    elif searchtype == "subclusters":
        blastoutputfile = open(outputfoldername + os.sep + "subclusterblastoutput.txt","w")
    elif searchtype == "knownclusters":
        blastoutputfile = open(outputfoldername + os.sep + "knownclusterblastoutput.txt","w")
    blastoutputfile.write(blastoutput)
    blastoutputfile.close()

class Score(object):
    __slots__ = ("hits", "core_gene_hits", "blast_score", "synteny_score",
                 "core_bonus", "scored_pairings")
    def __init__(self):
        self.hits = 0
        self.core_gene_hits = 0
        self.blast_score = 0
        self.synteny_score = 0
        self.core_bonus = 0
        self.scored_pairings = []

    def sort_score(self):
        """ the pre-existing, unexplained sort weighting """
        if self.core_gene_hits:
            self.core_bonus = 3
        return (self.hits
                + self.core_bonus
                + self.core_gene_hits / 100.
                + self.blast_score / 10e7
                + self.synteny_score)

def parse_clusterblast_dict(queries, clusters, cluster_name, allcoregenes):
    result = Score()
    hitpositions = []
    hitposcorelist = []
    cluster_locii = clusters[cluster_name][0]
    for query in queries:
        querynrhits = 0
        for subject in query.get_subjects_by_cluster(cluster_name):
            assert cluster_name == subject.genecluster
            if subject.locus_tag not in cluster_locii:
                continue
            index_pair = [query.index, cluster_locii.index(subject.locus_tag)]
            if index_pair in hitpositions:
                continue
            querynrhits += 1
            result.blast_score += subject.blastscore
            result.scored_pairings.append([query, subject])
            hitpositions.append(index_pair)
        if querynrhits:
            result.hits += 1
            hit_value = 0
            if query.id in allcoregenes:
                result.core_gene_hits += 1
                hit_value = 1
            hitposcorelist.extend([hit_value]*querynrhits)
    return result, hitpositions, hitposcorelist

def find_clusterblast_hitsgroups(hitpositions):
    #Find groups of hits
    hitgroupsdict = {}
    for p in hitpositions:
        if not hitgroupsdict.has_key(p[0]):
            hitgroupsdict[p[0]] = [p[1]]
        else:
            hitgroupsdict[p[0]].append(p[1])
    return hitgroupsdict

def calculate_synteny_score(hitgroupsdict, hitpositions, hitposcorelist):
    scored_queries = set()
    scored_hits = set()
    # Calculate synteny score; give score only if more than one hits
    # (otherwise no synteny possible), and only once for every query gene and every hit gene
    synteny_score = 0
    for n, p in enumerate(hitpositions[:-1]):
        n += 1 # we want a 1-index, not 0-index
        q = hitpositions[n]
        # Check if a gene homologous to this gene has already been
        # scored for synteny in the previous entry
        if p[1] in hitgroupsdict[q[0]] or p[0] in scored_queries or p[1] in scored_hits:
            continue
        if (abs(p[0] - q[0]) < 2) and abs(p[0]-q[0]) == abs(p[1]-q[1]):
            synteny_score += 1
            if hitposcorelist[n - 1] == 1 or hitposcorelist[n] == 1:
                synteny_score += 1
            scored_queries.add(p[0])
            scored_hits.add(p[1])
    return synteny_score

def score_clusterblast_output(clusters, allcoregenes, cluster_names_to_queries):
    #Score BLAST output on all gene clusters
    #Rank gene cluster hits based on 1) number of protein hits covering >25% sequence length or at least 100aa alignment, with >30% identity and 2) cumulative blast score
    #Find number of protein hits and cumulative blast score for each gene cluster
    results = {}
    for cluster_name, queries in cluster_names_to_queries.items():
        result, hitpositions, hitposcorelist = parse_clusterblast_dict(queries, clusters, cluster_name, allcoregenes)
        if result.hits <= 1:
            continue
        hitgroupsdict = find_clusterblast_hitsgroups(hitpositions)
        result.synteny_score = calculate_synteny_score(hitgroupsdict, hitpositions, hitposcorelist)
        # ensure at least two different subjects were found
        initial = hitpositions[0][1]
        for _, subject in hitpositions[1:]:
            if subject != initial:
                results[cluster_name] = result
                break
    #Sort gene clusters by score
    return sorted(results.items(), reverse=True, key=lambda x: x[1].sort_score())

def get_output_dir(options, searchtype):
    output_dir = None
    if searchtype == "general":
        options.clusterblast_outputfolder = options.full_outputfolder_path + os.sep + "clusterblast"
        output_dir = options.clusterblast_outputfolder
    elif searchtype == "subclusters":
        options.subclusterblast_outputfolder = options.full_outputfolder_path + os.sep + "subclusterblast"
        output_dir = options.subclusterblast_outputfolder
    elif searchtype == "knownclusters":
        options.knownclusterblast_outputfolder = options.full_outputfolder_path + os.sep + "knownclusterblast"
        output_dir = options.knownclusterblast_outputfolder
    else:
        raise ValueError("Invalid search type (%s) in clusterblast" % (searchtype))

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    return output_dir


def write_clusterblast_output(options, seq_record,clusterblastStorage, searchtype="general"):

    clusternumber = clusterblastStorage.clusternumber
    queryclusterprots = clusterblastStorage.queryclusterprots
    clusters = clusterblastStorage.clusters
    ranking = clusterblastStorage.ranking
    proteins = clusterblastStorage.proteins

    #Output for each hit: table of genes and locations of input cluster, table of genes and locations of hit cluster, table of hits between the clusters
    currentdir = os.getcwd()
    os.chdir(get_output_dir(options, searchtype))

    out_file = open("cluster" + str(clusternumber) + ".txt","w")
    out_file.write("ClusterBlast scores for " + seq_record.id + "\n")
    out_file.write("\nTable of genes, locations, strands and annotations of query cluster:\n")
    feature_by_id = utils.get_feature_dict_protein_id(seq_record)
    for i in queryclusterprots:
        cds = feature_by_id[i]
        if cds.strand == 1:
            strand = "+"
        else:
            strand = "-"
        out_file.write("\t".join([i, str(cds.location.nofuzzy_start), str(cds.location.nofuzzy_end), strand, utils.get_gene_annotation(cds)]) + "\t\n")
    out_file.write("\n\nSignificant hits: \n")
    top_hits = ranking[:100]
    for n, cluster_and_result in enumerate(top_hits):
        cluster = cluster_and_result[0]
        out_file.write("{}. {}\t{}\n".format(n + 1, cluster, clusters[cluster][1]))

    out_file.write("\n\nDetails:")
    for n, cluster_and_result in enumerate(top_hits):
        cluster, result = cluster_and_result
        # TODO: change to just result.hits during next minor version bump
        nrhits = result.hits + result.synteny_score + result.core_bonus
        out_file.write("\n\n>>\n")
        out_file.write("{}. {}\n".format(n + 1, cluster))
        out_file.write("Source: {}\n".format(clusters[cluster][1]))
        out_file.write("Type: {}\n".format(clusters[cluster][2]))
        out_file.write("Number of proteins with BLAST hits to this cluster: %d\n" % nrhits)
        out_file.write("Cumulative BLAST score: %d\n\n" % result.blast_score)
        out_file.write("Table of genes, locations, strands and annotations of subject cluster:\n")
        clusterproteins = clusters[cluster][0]
        for protein_name in clusterproteins:
            protein = proteins.get(protein_name)
            if protein:
                out_file.write(str(protein))
        out_file.write("\nTable of Blast hits (query gene, subject gene, %identity, blast score, %coverage, e-value):\n")
        if result.scored_pairings:
            for query, subject in result.scored_pairings:
                # TODO : check the trailing \t is meaningful
                out_file.write("{}\t{}\t\n".format(query.id, subject.get_table_string()))
        else:
            out_file.write("data not found\n")
        out_file.write("\n")
    out_file.close()
    os.chdir(currentdir)

def perform_clusterblast(options, seq_record, clusters, proteins):
    #Run BLAST on gene cluster proteins of each cluster and parse output
    geneclusters = utils.get_sorted_cluster_features(seq_record)
    debug_path = os.path.abspath(os.path.join(options.dbgclusterblast, "clusterblastoutput.txt"))
    with TemporaryDirectory(change=True) as tempdir:
        all_names, all_seqs, all_prots = [], [] ,[]
        prots_by_cluster = []
        for genecluster in geneclusters:
            names, seqs, prots = create_blast_inputs(genecluster, seq_record)
            all_names.extend(names)
            all_seqs.extend(seqs)
            all_prots.extend(prots)
            prots_by_cluster.append(prots)
        if options.dbgclusterblast and os.path.exists(debug_path):
            logging.debug("Skipping DIAMOND calculations, using results from %s instead",
                          debug_path)
            with open(debug_path, "r") as fh:
                blastoutput = fh.read()
            logging.debug("    Parsing results from given file...")
        else:
            logging.debug("Running DIAMOND gene cluster search..")
            utils.writefasta(all_names, all_seqs, "input.fasta")
            out, err, retcode = run_diamond("input.fasta", path.join(options.clusterblastdir, "geneclusterprots"), tempdir, options)
            if retcode != 0:
                logging.error("Running diamond failed: returned %s, stderr: %r, stdout: %r", retcode, err, out)
            logging.debug("   DIAMOND search finished. Parsing results...")

            with open("input.out", 'r') as fh:
                blastoutput = fh.read()

        write_raw_clusterblastoutput(options.full_outputfolder_path, blastoutput)


        minseqcoverage = 10
        minpercidentity = 30
        clusters_by_number, _ = parse_all_clusters(blastoutput, minseqcoverage,
                                                                 minpercidentity, seq_record)
        
        clusterblastStorage = utils.Storage()
        clusterblastStorage.clusters = clusters
        clusterblastStorage.proteins = proteins

        for genecluster, queryclusterprots in zip(geneclusters, prots_by_cluster):
            clusternumber = utils.get_cluster_number(genecluster)
            cluster_names_to_queries = clusters_by_number.get(clusternumber, {})
            allcoregenes = [utils.get_gene_acc(cds) for cds in utils.get_secmet_cds_features(seq_record)]
            ranking = score_clusterblast_output(clusters, allcoregenes, cluster_names_to_queries)

            # store all clusterblast related data in a utils.Storage object
            clusterblastStorage.clusternumber = clusternumber
            clusterblastStorage.queryclusterprots = queryclusterprots
            clusterblastStorage.ranking = ranking

            write_clusterblast_output(options, seq_record, clusterblastStorage)
