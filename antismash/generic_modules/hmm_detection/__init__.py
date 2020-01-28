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

"""General HMM detection module

"""
import logging
from os import path
from antismash.generic_modules import Signature
from antismash import config
from antismash import utils
from Bio.SeqFeature import SeqFeature, FeatureLocation

name = "hmmdetection"
short_description = name.capitalize()
priority = 5000



class HmmSignature(Signature):
    """HMM signature"""
    def __init__(self, name, description, cutoff, hmm_file):
        self.hmm_file = utils.get_full_path(__file__, hmm_file)
        self.name = name
        super(HmmSignature, self).__init__(name, 'model',
              description, cutoff, utils.get_full_path(__file__, hmm_file))

short_description = name.capitalize()

# The tuple is the name of the binary and whether it is an optional requirement
_required_binaries = [
    ('hmmsearch', False),
    ('hmmpress', False)
]

_markov_model = 'bgc_seeds.hmm'

_binary_extensions = [
    '.h3f',
    '.h3i',
    '.h3m',
    '.h3p',
]

#Define all profiles using details from hmmdetails.txt
hmmdetails = [line.split("\t") for line in open(utils.get_full_path(__file__, "hmmdetails.txt"),"r").read().split("\n") if line.count("\t") == 3]
_signature_profiles = [HmmSignature(details[0], details[1], int(details[2]), details[3]) for details in hmmdetails]

def check_prereqs():
    failure_messages = []
    for binary_name, optional in _required_binaries:
        if utils.locate_executable(binary_name) is None and not optional:
            failure_messages.append("Failed to locate executable for %r" %
                                    binary_name)

    hmm_files = []
    # Check if hmmdetails.txt is readable and well-formatted
    lineno = 1
    for line in open(utils.get_full_path(__file__, "hmmdetails.txt"),"r"):
        if line.count("\t") != 3:
            failure_messages.append("Failed to use HMM profile from line %s due to misformatting:\n %r" %
                                    (lineno, line))
            continue
        hmm_files.append(line.split('\t')[3].strip())
        lineno += 1

    #Check if cluster_rules.txt is readable and well-formatted
    lineno = 1
    for line in open(utils.get_full_path(__file__, "cluster_rules.txt"),"r"):
        if line.count("\t") != 3:
            failure_messages.append("Failed to use cluster rules from the line %s due to misformatting:\n %r" %
                                    (lineno, line))

        lineno += 1

    hmm = utils.get_full_path(__file__, _markov_model)
    if utils.locate_file(hmm) is None:
        # try to generate file from all specified profiles in hmmdetails
        try:
            with open(hmm, 'w') as all_hmms_handle:
                for hmm_file in hmm_files:
                    with open(utils.get_full_path(__file__, hmm_file), 'r') as handle:
                        all_hmms_handle.write(handle.read())
        except OSError:
            failure_messages.append('Failed to generate file {!r}'.format(hmm))

    for ext in _binary_extensions:
        binary = "{}{}".format(hmm, ext)
        if utils.locate_file(binary) is None:
            _, err, retcode = utils.run_hmmpress(hmm)
            if retcode != 0:
                failure_messages.append('Failed to hmmpress {!r}: {!r}'.format(hmm, err))
            break

    return failure_messages


def get_supported_cluster_types():
    "Get a list of all supported cluster types"
    clustertypes = [line.split("\t")[0] for line in open(utils.get_full_path(__file__, 'cluster_rules.txt'), "r")]
    # skip first line containing the header
    return clustertypes[1:]


def find_clusters(seq_record, rulesdict):
    #Functions that detects the gene clusters based on the identified core genes
    features = utils.get_cds_features(seq_record)
    clusters = []
    cfg = config.get_config()
    clusternr = cfg.next_clusternr

    for feature in features:
        within_cutoff = False
        if ('sec_met' not in feature.qualifiers) or (len([feat for feat in feature.qualifiers['sec_met'] if "Type: " in feat]) <= 0):
            continue
        feature_start = min(feature.location.start, feature.location.end)
        feature_end = max(feature.location.start, feature.location.end)
        feature_type = [feat for feat in feature.qualifiers['sec_met'] if "Type: " in feat][0].partition("Type: ")[2]
        if feature_type == "none":
            continue
        feature_cutoff = max([rulesdict[value][1] for value in feature_type.split("-")])
        feature_extension = max([rulesdict[value][2] for value in feature_type.split("-")])
        cluster = None

        if len(clusters) > 0:
            cluster = clusters[-1]
            cluster_end = cluster.location.end
            # Check cutoff
            cutoff = max(cluster.qualifiers['cutoff'][0], feature_cutoff)
            cutoff = max(cutoff, cluster.qualifiers['extension'][0] + feature_extension)
            within_cutoff = feature_start <= cluster_end + cutoff

        if not within_cutoff:
            if len(clusters) > 0:
                # Finalize the last extended cluster
                cluster = clusters[-1]
                cluster.location = FeatureLocation(max(0, cluster.location.start - cluster.qualifiers['extension'][0]), min(len(seq_record), cluster.location.end + cluster.qualifiers['extension'][0]))
            # Create new cluster
            new_cluster = SeqFeature(FeatureLocation(feature_start, feature_end), type="cluster")
            new_cluster.qualifiers['note'] = ["Cluster number: " + str(clusternr)]
            new_cluster.qualifiers['cutoff'] = [feature_cutoff]
            new_cluster.qualifiers['extension'] = [feature_extension]
            new_cluster.qualifiers['product'] = [feature_type]
            clusters.append(new_cluster)
            cluster = clusters[-1]
            clusternr += 1

        # Update cluster
        cluster.location = FeatureLocation(min(cluster.location.start, feature_start), max(cluster.location.end, feature_end))
        cluster.qualifiers['cutoff'] = [max(cluster.qualifiers['cutoff'][0], feature_cutoff)]
        cluster.qualifiers['extension']  = [max(cluster.qualifiers['extension'][0], feature_extension)]
        cluster.qualifiers['product'] =  ["-".join(list(set(cluster.qualifiers['product'][0].split('-')) | set(feature_type.split('-'))))]
        if "-" in cluster.qualifiers['product'][0]:
            cluster.qualifiers['product'] = ["-".join([ct for ct in cluster.qualifiers['product'][0].split('-') if ct != "other"])]

    if len(clusters) > 0:
        # Finalize the last extended cluster
        cluster = clusters[-1]
        cluster.location = FeatureLocation(max(0, cluster.location.start - cluster.qualifiers['extension'][0]), min(len(seq_record), cluster.location.end + cluster.qualifiers['extension'][0]))
    for cluster in clusters:
        #Add a note to specify whether a cluster lies on the contig/scaffold edge or not
        if cluster.location.start == 0 or cluster.location.end == len(seq_record):
            cluster.qualifiers['contig_edge'] = "True"
        else:
            cluster.qualifiers['contig_edge'] = "False"

    seq_record.features.extend(clusters)
    cfg.next_clusternr = clusternr


def filter_results(results, results_by_id):
    #Filter results by comparing scores of different models (for PKS systems)
    for line in open(utils.get_full_path(__file__, "filterhmmdetails.txt"),"r").read().split("\n"):
        filterhmms = line.split(",")
        for cds in results_by_id.keys():
            cdsresults = results_by_id[cds]
            hmmhits = [hit.query_id for hit in cdsresults]
            #Check if multiple competing HMM hits are present
            competing_hits = set(hmmhits) & set(filterhmms)
            if len(competing_hits) > 1:
                #Identify overlapping hits
                overlapping_groups = []
                for hit in cdsresults:
                    for otherhit in [cdsresult for cdsresult in cdsresults if hit != cdsresult]:
                        overlap = len(set(range(hit.hit_start, hit.hit_end)) & set(range(otherhit.hit_start, otherhit.hit_end)))
                        if overlap > 20:
                            added = "n"
                            for group in overlapping_groups:
                                if hit in group and otherhit in group:
                                    added = "y"
                                    break
                                elif hit in group and otherhit not in group:
                                    group.append(otherhit)
                                    added = "y"
                                    break
                                elif hit not in group and otherhit in group:
                                    group.append(hit)
                                    added = "y"
                                    break
                            if added == "n":
                                overlapping_groups.append([hit, otherhit])
                #Remove worst-scoring of overlapping hits
                for group in overlapping_groups:
                    highestscore = max([hit.bitscore for hit in group])
                    hit_with_highestscore = group[[hit.bitscore for hit in group].index(highestscore)]
                    to_delete = [hit for hit in group if hit != hit_with_highestscore]
                    for res in [res for res in results]:
                        if res in to_delete:
                            del results[results.index(res)]
                            del results_by_id[cds][results_by_id[cds].index(res)]
    return results, results_by_id


def filter_result_multiple(results, results_by_id):
    #Filter multiple results of the same model within a gene
    for cds in results_by_id.keys():
        best_hit_scores = {}
        to_delete = []
        for hit in results_by_id[cds]:
            if (hit.query_id not in best_hit_scores) or (best_hit_scores[hit.query_id] < hit.bitscore):
                best_hit_scores[hit.query_id] = hit.bitscore
        for hit in results_by_id[cds]:
            if (hit.bitscore < best_hit_scores[hit.query_id]):
                to_delete.append(hit)
        for hit in to_delete:
            del results[results.index(hit)]
            del results_by_id[cds][results_by_id[cds].index(hit)]
            if len(results_by_id[cds]) < 1:
                del results_by_id[cds]
    return results, results_by_id


def filter_result_overlapping_genes(results, results_by_id, overlaps, feature_by_id):
    # filter results of overlapping genes (only gene with the best score can retain its result)
    filterhmm_list = []
    overlap_id_with_result = {}
    for line in open(utils.get_full_path(__file__, "filterhmmdetails.txt"),"r").read().split("\n"):
        filterhmms = line.split(",")
        if filterhmms not in filterhmm_list:
            filterhmm_list.append(filterhmms)
    for cds in results_by_id.keys():
        if overlaps[1][cds] not in overlap_id_with_result.keys():
            overlap_id_with_result[overlaps[1][cds]] = [cds]
        elif cds not in overlap_id_with_result[overlaps[1][cds]]:
            overlap_id_with_result[overlaps[1][cds]].append(cds)
    for overlap_id in overlap_id_with_result.keys():
        best_hit_scores = {}
        for cds in overlap_id_with_result[overlap_id]:
            for hit in results_by_id[cds]:
                feature = feature_by_id[hit.hit_id]
                if (hit.query_id not in best_hit_scores) or (best_hit_scores[hit.query_id] < abs(feature.location.end - feature.location.start)):
                    best_hit_scores[hit.query_id] = abs(feature.location.end - feature.location.start)
        for cds in overlap_id_with_result[overlap_id]:
            to_delete = []
            for hit in results_by_id[cds]:
                feature = feature_by_id[hit.hit_id]
                if (abs(feature.location.end - feature.location.start) < best_hit_scores[hit.query_id]):
                    to_delete.append(hit)
                else: # filter for filterhmmdetails.txt
                    for filterhmms in filterhmm_list:
                        if hit.query_id not in filterhmms:
                            continue
                        for similar_hit in filterhmms:
                            if similar_hit not in best_hit_scores.keys():
                                continue
                            if (abs(feature.location.end - feature.location.start) < best_hit_scores[similar_hit]):
                                to_delete.append(hit)
                                break
            for hit in to_delete:
                del results[results.index(hit)]
                del results_by_id[cds][results_by_id[cds].index(hit)]
                if len(results_by_id[cds]) < 1:
                    del results_by_id[cds]
    return results, results_by_id


def create_rules_dict(enabled_clustertypes):
    "Create a cluster rules dictionary from the cluster rules file"
    rulesdict = {}
    first = True
    #TODO: We should move all user-customizable files into config subdirectory; the rulefiles are redundant also in hmm_detection_dblookup
    for line in open(utils.get_full_path(__file__, "cluster_rules.txt"),"r"):
        # skip the first line with the legend
        if first:
            first = False
            continue
        parts = line.split('\t')
        if len(parts) < 3:
            continue
        key = parts.pop(0)
        if key not in enabled_clustertypes:
            continue
        rules = parts.pop(0)
        cutoff = int(parts.pop(0)) * 1000
        extension = int(parts.pop(0)) * 1000
        rulesdict[key] = (rules, cutoff, extension)
    return rulesdict

def apply_cluster_rules(results_by_id, feature_by_id, enabled_clustertypes, rulesdict, overlaps):
    "Apply cluster rules to determine if HMMs lead to secondary metabolite core gene detection"
    typedict = {}
    cds_with_hits = sorted(results_by_id.keys(), key = lambda gene_id: feature_by_id[gene_id].location.start)
    for cds in cds_with_hits:
        _type = "none"
        #if typedict[cds] exist (the case of in-advance assignment from neighboring genes), use that instead of "none"
        if cds in typedict.keys():
            _type = typedict[cds]
        cdsresults = [res.query_id for res in results_by_id[cds]]

        for clustertype in enabled_clustertypes:
            if clustertype in _type:
                continue
            single_rules = [rule for rule in rulesdict[clustertype][0].split(" or ") if " & " not in rule and "cluster(" not in rule]
            combined_rules = [rule[1:-1] for rule in rulesdict[clustertype][0].split(" or ") if " & " in rule]
            cluster_rules = [rule[8:-1] for rule in rulesdict[clustertype][0].split(" or ") if "cluster(" in rule]
            minimum_rules = [rule[8:-1] for rule in rulesdict[clustertype][0].split(" or ") if "minimum(" in rule]
            if "-" in clustertype:
                cutoff = max([rulesdict[value][1] for value in clustertype.split("-")])
            else:
                cutoff = rulesdict[clustertype][1]
            #Assign cluster type if a single argument rule matches
            #Example rule format: "Domain1"
            if len(set(single_rules) & set(cdsresults)) >= 1:
                if not (_type != "none" and clustertype == "other"):
                    if _type == "none" or _type == "other" or _type == clustertype:
                        _type = clustertype
                    elif clustertype not in _type:
                        _type = clustertype + "-" + _type
                if _type != "other":
                    continue
            #Assign cluster type if a combinatorial argument rule matches
            #Example rule format: "(Domain1 & Domain2)"
            for rule in combined_rules:
                required_matches = rule.split(" & ")
                if len(set(required_matches) & set(cdsresults)) == len(required_matches):
                    if not (_type != "none" and clustertype == "other"):
                        if _type == "none" or _type == "other" or _type == clustertype:
                            _type = clustertype
                        elif clustertype not in _type:
                            _type = clustertype + "-" + _type
            if _type == clustertype and _type != "other":
                continue
            #Assign cluster type if distance-based combinatorial parameter matches
            #Example rule format: "cluster(Domain1,Domain2)"
            for rule in cluster_rules:
                #Find which needed domains are already found in present CDS
                required_matches = rule.split(",")
                cluster_results = list(set(required_matches) & set(cdsresults))
                missing_results = list(set(required_matches) - set(cdsresults))
                #If more than one, search nearby CDSs for the other needed domains
                if len(cluster_results) > 0:
                    locations = [feature_by_id[cds].location.start, feature_by_id[cds].location.end]
                    for othercds in results_by_id.keys():
                        needed_hits_found = list(set([res.query_id for res in results_by_id[othercds]]) & set(missing_results))
                        if len(needed_hits_found) > 0:
                            feature = feature_by_id[othercds]
                            flocations = [feature.location.start, feature.location.end]
                            #If hit found in nearby CDS, add it to the set of found domains relevant to the present rule
                            if min([abs(max(locations) - min(flocations)), abs(max(locations) - min(flocations))]) < cutoff:
                                cluster_results.extend(needed_hits_found)
                #If the rule is completely forfilled (all domains have a match in the same neighbourhood), assign cluster type
                if len(list(set(required_matches) - set(cluster_results))) == 0:
                    if not (_type != "none" and clustertype == "other"):
                        if _type == "none" or _type == "other" or _type == clustertype:
                            _type = clustertype
                        elif clustertype not in _type:
                            _type = clustertype + "-" + _type
                        break
            #Assign cluster type if distance-based combinatorial parameter matches with minimum number of domain hits
            #Example rule format: "minimum(4,[Domain1,Domain2], [Domain1])"
            #This denotes that at least 4 hits of either Domain1 or Domain2 should be found in the same region with at least 1 hit from Domain 1
            for rule in minimum_rules:
                #Find which needed domains are already found in present CDS
                min_number = int(rule.partition(",")[0])
                required_matches = rule.partition("[")[2].partition("]")[0].split(",")
                essential_matches = rule.rpartition("[")[2].partition("]")[0].split(",")
                if essential_matches == ['']:
                    essential_matches = []
                cluster_results = list(set(required_matches) & set(cdsresults))
                missing_essential = list(set(essential_matches) - set(calc_essential_found(essential_matches, cdsresults)))
                neighborcds = [] #neighboring cds that together met the requirements
                nrcds = 0
                if len(cluster_results) > 0:
                    nrcds = 1
                #If essentials found but havent fulfilled, or required found but havent fulfilled, search nearby CDSs for the other needed domains
                locations = [feature_by_id[cds].location.start, feature_by_id[cds].location.end]
                if ((len(list(set(essential_matches) & set(missing_essential))) > 0) or (nrcds > 0 and nrcds < min_number)):
                    #scan neighboring genes forward
                    for othercds in cds_with_hits:
                        if (overlaps[1][othercds] <= overlaps[1][cds]) or (overlaps[1][othercds] in [overlaps[1][ncds] for ncds in neighborcds]):
                            continue
                        feature = feature_by_id[othercds]
                        flocations = [feature.location.start, feature.location.end]
                        cluster_start = min(locations)
                        cluster_end = max(locations)
                        if ((min(flocations) >= cluster_start and max(flocations) <= cluster_end)
                                or (cluster_start - cutoff <= max(flocations) <= cluster_start)
                                or (cluster_end <= min(flocations) <= cluster_end + cutoff)):
                            othercds_results = [res.query_id for res in results_by_id[othercds]]
                            #check essential hits
                            essential_found = calc_essential_found(missing_essential, othercds_results)
                            if len(essential_found) > 0:
                                missing_essential = list(set(missing_essential) - set(essential_found))
                                #append to neighborcds for in-advance tagging
                                if not (othercds in neighborcds):
                                    neighborcds.append(othercds)
                                    locations.extend(flocations)
                            #check required hits
                            needed_hits_found = list(set(othercds_results) & set(required_matches))
                            if len(needed_hits_found) > 0:
                                nrcds += 1
                                #append to neighborcds for in-advance tagging
                                if not (othercds in neighborcds):
                                    neighborcds.append(othercds)
                                    locations.extend(flocations)
                        else:
                            break
                #If the rule is completely forfilled (all domains have a match in the same neighbourhood), assign cluster type
                #for the cds and all the corresponding neighbor cdss
                if nrcds >= min_number and len(missing_essential) == 0:
                    if not (_type != "none" and clustertype == "other"):
                        if _type == "none" or _type == "other" or _type == clustertype:
                            _type = clustertype
                        elif clustertype not in _type.split("-"):
                            _type = clustertype + "-" + _type
                    for ncds in neighborcds:
                        _ntype = "none"
                        if ncds in typedict.keys():
                            _ntype = typedict[ncds]
                        if not (_ntype != "none" and clustertype == "other"):
                            if _ntype == "none" or _ntype == "other" or _ntype == clustertype:
                                _ntype = clustertype
                            elif clustertype not in _ntype.split("-"):
                                _ntype = clustertype + "-" + _ntype
                        typedict[ncds] = _ntype
                    break
        #Save type to typedict
        typedict[cds] = _type
    return typedict


def detect_signature_genes(seq_record, enabled_clustertypes, options):
    "Function to be executed by module"
    feature_by_id = utils.get_feature_dict(seq_record)
    full_fasta = utils.get_multifasta(seq_record)
    rulesdict = create_rules_dict(enabled_clustertypes)
    results = []
    sig_by_name = {}
    results_by_id = {}
    for sig in _signature_profiles:
        sig_by_name[sig.name] = sig

    runresults = utils.run_hmmsearch(utils.get_full_path(__file__, 'bgc_seeds.hmm'), full_fasta, use_tempfile=True)
    for runresult in runresults:
        acc = runresult.accession.split('.')[0]
        # Store result if it is above cut-off
        for hsp in runresult.hsps:
            if hsp.query_id in sig_by_name:
                sig = sig_by_name[hsp.query_id]
            elif acc in sig_by_name:
                sig = sig_by_name[acc]
            else:
                logging.error('BUG: Failed to find signature for ID %s / ACC %s', hsp.query_id, acc)
                continue
            if hsp.bitscore > sig.cutoff:
                results.append(hsp)
                if hsp.hit_id not in results_by_id:
                    results_by_id[hsp.hit_id] = [hsp]
                else:
                    results_by_id[hsp.hit_id].append(hsp)

    #Get overlap tables (for overlap filtering etc)
    overlaps = utils.get_overlaps_table(seq_record)

    #Filter results by comparing scores of different models (for PKS systems)
    results, results_by_id = filter_results(results, results_by_id)

    # Filter results of overlapping genes (only for plants)
    if options.taxon == 'plants':
        results, results_by_id = filter_result_overlapping_genes(results, results_by_id, overlaps, feature_by_id)

    #Filter multiple results of the same model in one gene
    results, results_by_id = filter_result_multiple(results, results_by_id)

    #Use rules to determine gene clusters
    typedict = apply_cluster_rules(results_by_id, feature_by_id, enabled_clustertypes, rulesdict, overlaps)

    #Find number of sequences on which each pHMM is based
    nseqdict = get_nseq()

    #Save final results to seq_record
    for cds in results_by_id.keys():
        feature = feature_by_id[cds]
        _update_sec_met_entry(feature, results_by_id[cds], typedict[cds], nseqdict)

    find_clusters(seq_record, rulesdict)

    #Find additional NRPS/PKS genes in gene clusters
    add_additional_nrpspks_genes(typedict, results_by_id, seq_record, nseqdict)
    #Add details of gene cluster detection to cluster features
    store_detection_details(results_by_id, rulesdict, seq_record)
    #If all-orfs option on, remove irrelevant short orfs
    if options.all_orfs:
        remove_irrelevant_allorfs(seq_record)

def get_nseq():
    nseqdict = {}
    for hmm in _signature_profiles:
        hmmfile = hmm.hmm_file
        for line in open(hmmfile, 'r'):
            if line.startswith('NSEQ '):
                nseqdict[hmm.name] = line[6:].strip()
                break
        if not hmm.name in nseqdict:
            nseqdict[hmm.name] = "?"

    return nseqdict

def overlaps(feature1, feature2):
    if (feature2.location.start <= feature1.location.start <= feature2.location.end) or (feature2.location.start <= feature1.location.end <= feature2.location.end):
        return True
    else:
        return False

def remove_irrelevant_allorfs(seq_record):
    #Get features
    allfeatures = utils.get_cds_features(seq_record)
    #Remove auto-orf features without unique sec_met qualifiers; remove glimmer ORFs overlapping with sec_met auto-orfs not catched by Glimmer
    auto_orf_features = [feature for feature in allfeatures if feature.qualifiers.has_key('note') and "auto-all-orf" in feature.qualifiers['note']]
    other_features = [feature for feature in allfeatures if not feature.qualifiers.has_key('note') or "auto-all-orf" not in feature.qualifiers['note']]
    to_delete = []
    for autofeature in auto_orf_features:
        if not autofeature.qualifiers.has_key("sec_met"):
            to_delete.append(autofeature)
        else:
            glimmer_has_sec_met = False
            for otherfeature in other_features:
                if overlaps(autofeature, otherfeature) and otherfeature.qualifiers.has_key('sec_met'):
                    to_delete.append(autofeature)
                    glimmer_has_sec_met = True
            if glimmer_has_sec_met == False:
                for otherfeature in other_features:
                    if overlaps(autofeature, otherfeature) and not otherfeature.qualifiers.has_key('sec_met'):
                        to_delete.append(otherfeature)
    featurenrs = []
    idx = 0
    for feature in seq_record.features:
        if feature in to_delete:
            featurenrs.append(idx)
        idx += 1
    featurenrs.reverse()
    for featurenr in featurenrs:
        del seq_record.features[featurenr]

def add_additional_nrpspks_genes(typedict, results_by_id, seq_record, nseqdict):
    nrpspksdomains = ["PKS_KS", "PKS_AT", "ATd", "ene_KS", "mod_KS", "hyb_KS", "itr_KS", "tra_KS", "Condensation", "AMP-binding", "A-OX"]
    clustercdsfeatures = utils.get_withincluster_cds_features(seq_record)
    othercds_with_results = [cds for cds in clustercdsfeatures if results_by_id.has_key(utils.get_gene_id(cds)) and typedict[utils.get_gene_id(cds)] == "none"]
    for cds in othercds_with_results:
        cdsresults = [res.query_id for res in results_by_id[utils.get_gene_id(cds)]]
        if len(set(nrpspksdomains) & set(cdsresults)) >= 1:
            _update_sec_met_entry(cds, results_by_id[utils.get_gene_id(cds)], "other", nseqdict)


def store_detection_details(results_by_id, rulesdict, seq_record):
    clusters = utils.get_cluster_features(seq_record)
    for cluster in clusters:
        type_combo = utils.get_cluster_type(cluster)
        if '-' in type_combo:
            clustertypes = type_combo.split('-')
        else:
            clustertypes = [type_combo]

        if not 'note' in cluster.qualifiers:
            cluster.qualifiers['note'] = []
        rule_string = "Detection rule(s) for this cluster type:"
        for clustertype in clustertypes:
            rule_string += " %s: (%s);" % (clustertype, rulesdict[clustertype][0])

        cluster.qualifiers['note'].append(rule_string)


def calc_essential_found(essential_matches, cdsresults):
    "Given array of essential matches, returns array of essential matches forfilled by the cdsresults"
    essential_found = []
    for match in essential_matches:
        for domain in match.split("/"):
            if domain in cdsresults:
                essential_found.append(match)
                break
    return essential_found


def _update_sec_met_entry(feature, results, clustertype, nseqdict):
    result = "; ".join(["%s (E-value: %s, bitscore: %s, seeds: %s)" % (res.query_id, res.evalue, res.bitscore, nseqdict.get(res.query_id, '?'))  for res in results])

    if not 'sec_met' in feature.qualifiers:
        feature.qualifiers['sec_met'] = [
            "Type: %s" % clustertype,
            "Domains detected: %s" % (result),
            "Kind: biosynthetic"
        ]
    else:
        for i, ann in enumerate(feature.qualifiers['sec_met']):
            if not ann.startswith('Domains detected: '):
                continue
            feature.qualifiers['sec_met'][i] = ann + "; %s" % (result)


__all__ = ["check_prereqs", "detect_signature_genes"]
