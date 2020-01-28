# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2016 Thomas Wolf
# Leibniz Institute for Natural Product Research and
# Infection Biology -- Hans-Knoell-Institute (HKI)
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Implementation of the CASSIS method for the motif-based prediction of SM gene clusters"""

import logging
import os
import csv
from xml.etree import cElementTree as ElementTree
import shutil

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
from Bio.SeqFeature import FeatureLocation

from antismash import utils


class AntiSmashError(RuntimeError):
    '''Base class for all antiSMASH custom errors'''
    pass


class InvalidLocationError(AntiSmashError):
    '''Thrown when running into invalid gene locations during runtime'''
    pass


class DuplicatePromoterError(AntiSmashError):
    '''Thrown when running into valid but duplicate promoter sequences during runtime'''
    pass


### plugin variables ###

name = "cassis"
short_description = name + ": Detect secondary metabolite gene cluster (motif based)"
priority = 2 # first run hmmdetect plugin to detect core genes (anchor genes) --> seed for cluster prediction with cassis

_required_binaries = [
    ("meme", "4.11.1"),
    ("fimo", "4.11.1"),
]

# all possible promoter sets for motif detection
# plus --> include <plus> promoters downstream the anchor gene's promoter
# minus --> include <minus> promoters upstream the anchor gene's promoter
_plus_minus = []
for plus in range(16):
    for minus in range(16):
        if plus + minus + 1 >= 4:
            _plus_minus.append({"plus": plus, "minus": minus})


### utility methods ###
def check_prereqs(options):
    """Check for prerequisites"""
    failure_messages = []
    for binary_name, binary_version in _required_binaries:
        if utils.locate_executable(binary_name) is None:
            failure_messages.append("Failed to locate executable for {!r}".format(binary_name))
        # TODO: Check binary version here

    return failure_messages


def get_versions(options):
    """Get all utility versions"""
    return map(lambda x: " ".join(x), _required_binaries)


### helper methods ###
def mprint(plus, minus):
    """Motif-print: nicely format motif name in plus/minus style"""
    return "+{:02d}_-{:02d}".format(plus, minus)


def get_promoter_id(promoter):
    """Return promoter ID string dependend on involved gene(s)"""
    if len(promoter["id"]) == 1: # 1 gene --> 1 promoter
        return promoter["id"][0]
    else: # 2 bidirectional genes --> 1 shared promoter
        return "{}+{}".format(promoter["id"][0], promoter["id"][1])


### main method ###
def detect(seq_record, options):
    """Use core genes (anchor genes) from hmmdetect as seeds to detect gene clusters"""
    logging.info("Detecting gene clusters using CASSIS")

    # TODO options? cassis settings/parameters?

    # get core genes from hmmdetect --> necessary CASSIS input, aka "anchor genes"
    anchor_genes = get_anchor_genes(seq_record)
    logging.debug("Record has {} anchor genes".format(len(anchor_genes)))
    if len(anchor_genes) == 0:
        return

    # filter all genes in record for neighbouring genes with overlapping annotations
    genes = utils.get_all_features_of_type(seq_record, "gene")
    logging.debug("Record has {} features of type 'gene'".format(len(genes)))
    if len(genes) == 0:
        return
    genes, ignored_genes = ignore_overlapping(genes)

    # compute promoter sequences/regions --> necessary for motif prediction (MEME and FIMO input)
    promoters = []
    try:
        # why these values? see "Wolf et al (2015): CASSIS and SMIPS ..."
        upstream_tss = 1000  # nucleotides upstream TSS
        downstream_tss = 50  # nucleotides downstream TSS
        promoters = get_promoters(seq_record, genes, upstream_tss, downstream_tss, options)
    except (InvalidLocationError, DuplicatePromoterError):
        logging.warning("CASSIS discovered an error while working on the promoter sequences, skipping CASSIS analysis")
        return

    if len(promoters) == 0:
        logging.debug("CASSIS found zero promoter regions, skipping CASSIS analysis")
        return
    if len(promoters) < 3:
        logging.debug("Sequence {!r} yields less than 3 promoter regions, skipping CASSIS analysis".format(seq_record.name))
        return
    if len(promoters) < 40:
        logging.debug("Sequence {!r} yields only {} promoter regions".format(seq_record.name, len(promoters)))
        logging.debug("Cluster detection on small sequences may lead to incomplete cluster predictions")

    store_promoters(promoters, seq_record)

    cluster_predictions = {} # {anchor gene : cluster predictions}
    for i in xrange(len(anchor_genes)):
        anchor = anchor_genes[i]
        logging.debug("Detecting cluster around anchor gene {!r} ({} of {})".format(anchor, i + 1, len(anchor_genes)))

        anchor_promoter = get_anchor_promoter(anchor, promoters)
        if anchor_promoter is None:
            logging.debug("No promoter region for {!r}, skipping this anchor gene".format(anchor))
            continue

        # predict motifs with MEME ("de novo")
        meme_dir = os.path.join(options.outputfoldername, "meme", anchor)
        promoter_sets = get_promoter_sets(meme_dir, anchor_promoter, promoters)
        exit_code = utils.run_meme(meme_dir, options)
        if exit_code != 0:
            logging.warning("MEME discovered a problem (exit code {}), skipping this anchor gene".format(exit_code))
            continue
        motifs = filter_meme_results(meme_dir, promoter_sets, anchor)

        if len(motifs) == 0:
            logging.debug("Could not predict motifs around {!r}, skipping this anchor gene".format(anchor))
            continue

        # search motifs with FIMO ("scanning")
        fimo_dir = os.path.join(options.outputfoldername, "fimo", anchor)
        exit_code = utils.run_fimo(meme_dir, fimo_dir, seq_record, options)
        if exit_code != 0:
            logging.warning("FIMO discovered a problem (exit code {}), skipping this anchor gene".format(exit_code))
            continue
        motifs = filter_fimo_results(motifs, fimo_dir, promoters, anchor_promoter, options)

        if len(motifs) == 0:
            logging.debug("Could not find motif occurrences for {!r}, skipping this anchor gene".format(anchor))
            continue

        # TODO SiTaR (http://bioinformatics.oxfordjournals.org/content/27/20/2806):
        # Alternative to MEME and FIMO. Part of the original CASSIS implementation.
        # No motif prediction (no MEME). Motif search with SiTaR (instead if FIMO).
        # Have to provide a file in FASTA format with binding site sequences of at least one transcription factor.
        # Will result in binding sites per promoter (like FIMO) --> find islands
        #
        # implement: YES? NO?

        # find islands of binding sites around anchor gene
        islands = get_islands(anchor_promoter, motifs, promoters, options)
        logging.debug("{} cluster predictions for {!r}".format(len(islands), anchor))

        # return cluster predictions sorted by border abundance
        # (most abundant --> "best" prediction)
        cluster_predictions[anchor] = sort_by_abundance(islands)
        cluster_predictions[anchor] = check_cluster_predictions(cluster_predictions[anchor], seq_record, promoters, ignored_genes)

        store_clusters(anchor, cluster_predictions[anchor], seq_record)

    if not options.skip_cleanup:
        logging.debug("Cleaning up MEME and FIMO output directories")
        cleanup_outdir(anchor_genes, cluster_predictions, options)


### additional methods ###
def get_anchor_genes(seq_record):
    """Return all genes which are putative cluster anchor genes"""
    anchor_genes = []

    for feature in seq_record.features:
        if "sec_met" in feature.qualifiers:
            for entry in feature.qualifiers['sec_met']:
                if entry.startswith('Type: ') and not entry.endswith('none'):
                    anchor_genes.append(utils.get_gene_id(feature))

    return anchor_genes


def ignore_overlapping(genes):
    """Ignore genes with overlapping locations (skip the second gene of an overlapping couple)"""
    ignored = []

    overlap = True
    while overlap: # check again until we didn't find any overlap in the entire (remaining) gene list
        overlap = False
        non_overlapping = [genes[0]]

        for i in xrange(1, len(genes)):
            if utils.features_overlap(genes[i-1], genes[i]):
                logging.debug("Ignoring {!r} (overlapping with {!r})".format(
                    utils.get_gene_id(genes[i]), utils.get_gene_id(genes[i-1])))
                ignored.append(genes[i])
                overlap = True
            else:
                non_overlapping.append(genes[i])

        genes = non_overlapping

    if ignored:
        logging.debug("Ignoring {} genes due to overlapping locations".format(len(ignored)))

    return (genes, ignored)


def get_promoters(seq_record, genes, upstream_tss, downstream_tss, options):
    """Compute promoter sequences for each gene in the sequence record"""
    logging.debug("Computing promoter sequences")

    min_promoter_length = 6
    max_promoter_length = (upstream_tss + downstream_tss) * 2 + 1

    record_seq_length = len(seq_record.seq)
    promoters = []
    invalid = 0

    pos_handle = open(os.path.join(options.outputfoldername, seq_record.name + "_promoter_positions.csv"), "w")
    pos_handle.write("\t".join(["#", "promoter", "start", "end", "length"]) + "\n")
    seq_handle = open(os.path.join(options.outputfoldername, seq_record.name + "_promoter_sequences.fasta"), "w")

    skip = 0  # helper var for shared promoter of bidirectional genes
    for i in xrange(len(genes)):

        if skip:  # two genes share the same promotor --> did computation with first gene, skip second gene
            skip = 0

        elif len(genes) == 1:  # only one gene within record

            if genes[i].location.strand == 1:
                # 1 (for explanation of these numbers see file promoterregions.png)
                if (genes[i].location.start - upstream_tss >= 0
                        and genes[i].location.end > genes[i].location.start + downstream_tss):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start - upstream_tss,
                        "end": genes[i].location.start + downstream_tss
                        # fuzzy (>|<) gene locations will be transformed to exact promoter locations
                        # we could save the fuzzy locations for promoters, too, via a FeatureLocation object
                        # but we use/calculate with the exact promoter locations anyway, here and later on
                    })
                # 2
                elif (genes[i].location.start - upstream_tss < 0
                        and genes[i].location.end > genes[i].location.start + downstream_tss):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": 0,
                        "end": genes[i].location.start + downstream_tss
                    })
                # 3
                elif (genes[i].location.start - upstream_tss >= 0
                        and genes[i].location.start + downstream_tss >= genes[i].location.end):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start - upstream_tss,
                        "end": genes[i].location.end
                    })
                # 7
                elif (genes[i].location.start - upstream_tss < 0
                        and genes[i].location.start + downstream_tss >= genes[i].location.end):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": 0,
                        "end": genes[i].location.end
                    })
                else:
                    logging.error("BUG: Problem with promoter of gene {!r}".format(utils.get_gene_id(genes[i])))
                    raise InvalidLocationError

            elif genes[i].location.strand == -1:
                # 4
                if (genes[i].location.start < genes[i].location.end - downstream_tss
                        and genes[i].location.end + upstream_tss <= record_seq_length):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": genes[i].location.end + upstream_tss
                    })
                # 5
                elif (genes[i].location.start < genes[i].location.end - downstream_tss
                        and genes[i].location.end + upstream_tss > record_seq_length):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": record_seq_length
                    })
                # 6
                elif (genes[i].location.start >= genes[i].location.end - downstream_tss
                        and genes[i].location.end + upstream_tss <= record_seq_length):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": genes[i].location.end + upstream_tss
                    })
                # 8
                elif (genes[i+1].location.start >= genes[i].location.end - upstream_tss
                        and genes[i].location.end + upstream_tss > record_seq_length):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": record_seq_length
                    })
                else:
                    logging.error("BUG: Problem with promoter of gene {!r}".format(utils.get_gene_id(genes[i])))
                    raise InvalidLocationError

        # first gene of the record AND NOT special case #9
        elif (i == 0 and not (genes[i].location.strand == -1
                                and genes[i+1].location.strand == 1
                                and genes[i].location.end + upstream_tss >= genes[i+1].location.start - upstream_tss)):

            if genes[i].location.strand == 1:
                #1
                if (genes[i].location.start - upstream_tss >= 0
                        and genes[i].location.end > genes[i].location.start + downstream_tss):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start - upstream_tss,
                        "end": genes[i].location.start + downstream_tss
                    })
                #2
                elif (genes[i].location.start - upstream_tss < 0
                        and genes[i].location.end > genes[i].location.start + downstream_tss):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": 0,
                        "end": genes[i].location.start + downstream_tss
                    })
                #3
                elif (genes[i].location.start - upstream_tss >= 0
                        and genes[i].location.start + downstream_tss >= genes[i].location.end):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start - upstream_tss,
                        "end": genes[i].location.end
                    })
                #7
                elif (genes[i].location.start - upstream_tss < 0
                        and genes[i].location.start + downstream_tss >= genes[i].location.end):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": 0,
                        "end": genes[i].location.end
                    })
                else:
                    logging.error("BUG: Problem with promoter of gene {!r}".format(utils.get_gene_id(genes[i])))
                    raise InvalidLocationError

            elif genes[i].location.strand == -1:
                #4
                if (genes[i].location.start < genes[i].location.end - downstream_tss
                        and genes[i+1].location.start > genes[i].location.end + upstream_tss):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": genes[i].location.end + upstream_tss
                    })
                #5
                elif (genes[i].location.start < genes[i].location.end - downstream_tss
                        and genes[i+1].location.start <= genes[i].location.end + upstream_tss):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": genes[i+1].location.start - 1
                    })
                #6
                elif (genes[i].location.start >= genes[i].location.end - downstream_tss
                        and genes[i+1].location.start > genes[i].location.end + upstream_tss):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": genes[i].location.end + upstream_tss
                    })
                #8
                elif (genes[i+1].location.start <= genes[i].location.end + upstream_tss
                        and genes[i].location.start >= genes[i].location.end - downstream_tss):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": genes[i+1].location.start - 1
                    })
                else:
                    logging.error("BUG: Problem with promoter of gene {!r}".format(utils.get_gene_id(genes[i])))
                    raise InvalidLocationError

        # last gene of record
        elif i == len(genes) - 1 and not skip:

            if genes[i].location.strand == 1:
                #1
                if (genes[i-1].location.end < genes[i].location.start - upstream_tss
                        and genes[i].location.end > genes[i].location.start + downstream_tss):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start - upstream_tss,
                        "end": genes[i].location.start + downstream_tss
                    })
                #2
                elif (genes[i-1].location.end >= genes[i].location.start - upstream_tss
                        and genes[i].location.end > genes[i].location.start + downstream_tss):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i-1].location.end + 1,
                        "end": genes[i].location.start + downstream_tss
                    })
                #3
                elif (genes[i-1].location.end < genes[i].location.start - upstream_tss
                        and genes[i].location.start + downstream_tss >= genes[i].location.end):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start - upstream_tss,
                        "end": genes[i].location.end
                    })
                #7
                elif (genes[i-1].location.end >= genes[i].location.start - upstream_tss
                        and genes[i].location.start + downstream_tss >= genes[i].location.end):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i-1].location.end + 1,
                        "end": genes[i].location.end
                    })
                else:
                    logging.error("BUG: Problem with promoter of gene {!r}".format(utils.get_gene_id(genes[i])))
                    raise InvalidLocationError

            elif genes[i].location.strand == -1:
                #4
                if (genes[i].location.start < genes[i].location.end - downstream_tss
                        and genes[i].location.end + upstream_tss <= record_seq_length):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": genes[i].location.end + upstream_tss
                    })
                #5
                elif (genes[i].location.start < genes[i].location.end - downstream_tss
                        and genes[i].location.end + upstream_tss > record_seq_length):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": record_seq_length
                    })
                #6
                elif (genes[i].location.start >= genes[i].location.end - downstream_tss
                        and genes[i].location.end + upstream_tss <= record_seq_length):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": genes[i].location.end + upstream_tss
                    })
                #8
                elif (genes[i+1].location.start <= genes[i].location.end + upstream_tss
                        and genes[i].location.end + upstream_tss > record_seq_length):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": record_seq_length
                    })
                else:
                    logging.error("BUG: Problem with promoter of gene {!r}".format(utils.get_gene_id(genes[i])))
                    raise InvalidLocationError

        # special-case 9
        elif (genes[i].location.strand == -1
                and genes[i+1].location.strand == 1
                and genes[i].location.end + upstream_tss >= genes[i+1].location.start - upstream_tss):

            #9 (1+4)
            if (genes[i].location.end > genes[i].location.start + downstream_tss
                    and genes[i].location.start < genes[i].location.end - downstream_tss):
                promoters.append({
                    "id": [utils.get_gene_id(genes[i]), utils.get_gene_id(genes[i+1])],
                    "start": genes[i].location.end - downstream_tss,
                    "end": genes[i+1].location.start + downstream_tss
                })
            #9 (3+4)
            elif (genes[i].location.start < genes[i].location.end - downstream_tss
                    and genes[i+1].location.start + downstream_tss >= genes[i+1].location.end):
                promoters.append({
                    "id": [utils.get_gene_id(genes[i]), utils.get_gene_id(genes[i+1])],
                    "start": genes[i].location.end - downstream_tss,
                    "end": genes[i+1].location.end
                })
            #9 (1+6)
            elif (genes[i].location.start >= genes[i].location.end - downstream_tss
                    and genes[i+1].location.end > genes[i+1].location.start + downstream_tss):
                promoters.append({
                    "id": [utils.get_gene_id(genes[i]), utils.get_gene_id(genes[i+1])],
                    "start": genes[i].location.start,
                    "end": genes[i+1].location.start + downstream_tss
                })
            #9 (3+6)
            elif (genes[i].location.start >= genes[i].location.end - downstream_tss
                    and genes[i+1].location.start + downstream_tss >= genes[i+1].location.end):
                promoters.append({
                    "id": [utils.get_gene_id(genes[i]), utils.get_gene_id(genes[i+1])],
                    "start": genes[i].location.start,
                    "end": genes[i+1].location.end
                })
            else:
                logging.error("BUG: Problem with promoter of gene {!r}".format(utils.get_gene_id(genes[i])))
                raise InvalidLocationError

            skip = 1

        # "normal" cases
        elif not skip:

            if genes[i].location.strand == 1:
                #1
                if (genes[i-1].location.end < genes[i].location.start - upstream_tss
                        and genes[i].location.end > genes[i].location.start + downstream_tss):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start - upstream_tss,
                        "end": genes[i].location.start + downstream_tss
                    })
                #2
                elif (genes[i-1].location.end >= genes[i].location.start - upstream_tss
                        and genes[i].location.end > genes[i].location.start + downstream_tss):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i-1].location.end + 1,
                        "end": genes[i].location.start + downstream_tss
                    })
                #3
                elif (genes[i-1].location.end < genes[i].location.start - upstream_tss
                        and genes[i].location.start + downstream_tss >= genes[i].location.end):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start - upstream_tss,
                        "end": genes[i].location.end
                    })
                #7
                elif (genes[i-1].location.end >= genes[i].location.start - upstream_tss
                        and genes[i].location.start + downstream_tss >= genes[i].location.end):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i-1].location.end + 1,
                        "end": genes[i].location.end
                    })
                else:
                    logging.error("BUG: Problem with promoter of gene {!r}".format(utils.get_gene_id(genes[i])))
                    raise InvalidLocationError

            elif genes[i].location.strand == -1:
                #4
                if (genes[i].location.start < genes[i].location.end - downstream_tss
                        and genes[i+1].location.start > genes[i].location.end + upstream_tss):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": genes[i].location.end + upstream_tss
                    })
                #5
                elif (genes[i].location.start < genes[i].location.end - downstream_tss
                        and genes[i+1].location.start <= genes[i].location.end + upstream_tss):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.end - downstream_tss,
                        "end": genes[i+1].location.start - 1
                    })
                #6
                elif (genes[i].location.start >= genes[i].location.end - downstream_tss
                        and genes[i+1].location.start > genes[i].location.end + upstream_tss):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": genes[i].location.end + upstream_tss
                    })
                #8
                elif (genes[i+1].location.start <= genes[i].location.end + upstream_tss
                        and genes[i].location.start >= genes[i].location.end - downstream_tss):
                    promoters.append({
                        "id": [utils.get_gene_id(genes[i])],
                        "start": genes[i].location.start,
                        "end": genes[i+1].location.start - 1
                    })
                else:
                    logging.error("BUG: Problem with promoter of gene {!r}".format(utils.get_gene_id(genes[i])))
                    raise InvalidLocationError

        # negative start position or stop position "beyond" record --> might happen in very small records
        if promoters[-1]["start"] < 0:
            promoters[-1]["start"] = 0
        if promoters[-1]["end"] > record_seq_length - 1:
            promoters[-1]["end"] = record_seq_length - 1

        # write promoter positions and sequences to file
        if not skip:
            promoter_sequence = seq_record.seq[promoters[-1]["start"]:promoters[-1]["end"]+1]
            promoter_length = len(promoter_sequence)

            invalid_promoter_sequence = ""

            # check if promoter length is valid
            if promoter_length < min_promoter_length or promoter_length > max_promoter_length:
                invalid_promoter_sequence = "length"

            # check if a, c, g and t occur at least once in the promoter sequence
            elif "A" not in promoter_sequence.upper():
                invalid_promoter_sequence = "A"
            elif "C" not in promoter_sequence.upper():
                invalid_promoter_sequence = "C"
            elif "G" not in promoter_sequence.upper():
                invalid_promoter_sequence = "G"
            elif "T" not in promoter_sequence.upper():
                invalid_promoter_sequence = "T"

            if invalid_promoter_sequence:
                invalid += 1

                if invalid_promoter_sequence == "length":
                    logging.warning("Promoter {!r} is invalid (length is {})".format(
                        get_promoter_id(promoters[-1]), promoter_length))
                else:
                    # especially SiTaR doesn't like such missings
                    logging.warning("Promoter {!r} is invalid (sequence without {!r})".format(
                        get_promoter_id(promoters[-1]), invalid_promoter_sequence))

                # more details for debug logging
                logging.debug("Invalid promoter {!r}\n start {}\n end {}\n length {}\n".format(
                    promoters[-1]["id"], promoters[-1]["start"], promoters[-1]["end"], promoter_length))

                promoters.pop() # remove last (invalid!) promoter

            else:
                promoters[-1]["seq"] = promoter_sequence

                # write promoter positions to file
                pos_handle.write("\t".join(map(str,
                    [len(promoters), get_promoter_id(promoters[-1]), promoters[-1]["start"] + 1, promoters[-1]["end"] + 1, promoter_length])) + "\n")

                # write promoter sequences to file
                SeqIO.write(
                    SeqRecord(promoter_sequence, id = get_promoter_id(promoters[-1]),
                    description = "length={}bp".format(promoter_length)),
                    seq_handle,
                    "fasta"
                )

        # check if promoter IDs are unique
        if len(promoters) >= 2 and get_promoter_id(promoters[-1]) == get_promoter_id(promoters[-2]):
            logging.error("Promoter {!r} occurs at least twice. This may be caused by overlapping gene annotations".format(
                get_promoter_id(promoters[-1])))
            raise DuplicatePromoterError

    if invalid:
        logging.debug("Ignoring {} promoters due to invalid promoter sequences".format(invalid))

    logging.debug("Found {} promoter sequences for {} genes".format(len(promoters), len(genes)))
    return promoters


def get_anchor_promoter(anchor, promoters):
    """Find the name of the promoter which includes the anchor gene"""
    # the promoter ID is not (necessarily) equal to the anchor ID!
    for i in xrange(len(promoters)):
        if anchor in promoters[i]["id"]:
            return i

    return None


def get_promoter_sets(meme_dir, anchor_promoter, promoters):
    """Prepare sets of promoter sequences and motif subdirectories"""
    promoter_sets = []

    if not os.path.exists(meme_dir):
        os.makedirs(meme_dir)

    # prepare sets of promoter sequences (MEME input)
    indices = set() # helper variable to monitor unique start_index/end_index tupels
    for pm in _plus_minus:
        start_index = anchor_promoter - pm["minus"]
        end_index = anchor_promoter + pm["plus"]

        if start_index < 0: # anchor promoter near beginning of record --> truncate
            logging.debug("Promoter set " + mprint(pm["plus"], pm["minus"]) + " exceeds upstream record border")
            start_index = 0

        if end_index > len(promoters) - 1: # anchor promoter near end of record --> truncate
            logging.debug("Promoter set " + mprint(pm["plus"], pm["minus"]) + " exceeds downstream record border")
            end_index = len(promoters) - 1

        # discard promoter sets, which reappear due to truncation
        if (start_index, end_index) not in indices:
            indices.add((start_index, end_index))

            # check (again, compare init of _plus_minus) if the promoter set has at least 4 promoters
            if end_index - start_index + 1 >= 4:
                promoter_sets.append({"plus": pm["plus"], "minus": pm["minus"], "score": None})

                pm_dir = os.path.join(meme_dir, mprint(pm["plus"], pm["minus"]))
                if not os.path.exists(pm_dir):
                    os.makedirs(pm_dir)

                # write promoter sequences to fasta file, in respective "plus-minus" subdir
                with open(os.path.join(pm_dir, "promoters.fasta"), "w") as pm_handle:
                    for i in xrange(start_index, end_index + 1):
                        seq = SeqRecord(promoters[i]["seq"],
                            id = get_promoter_id(promoters[i]),
                            description = "length={}bp".format(len(promoters[i]["seq"])))
                        if i == anchor_promoter: # mark anchor gene
                            seq.id += "__ANCHOR" # must be part of id, otherwise MEME woun't recognize it
                        SeqIO.write(seq, pm_handle, "fasta")
            else:
                logging.debug("Too short promoter set " + mprint(pm["plus"], pm["minus"]))
        else:
            logging.debug("Duplicate promoter set " + mprint(pm["plus"], pm["minus"]))

    return promoter_sets


def filter_meme_results(meme_dir, promoter_sets, anchor):
    """Analyse and filter MEME results"""
    for motif in promoter_sets:
        xml_file = os.path.join(meme_dir, mprint(motif["plus"], motif["minus"]), "meme.xml")
        e = ElementTree.parse(xml_file).getroot()
        reason = e.find("model/reason_for_stopping").text
        anchor_seq_id = ""

        # no motif found for given e-value cutoff :-(
        if "Stopped because motif E-value > " in reason:
            logging.debug("MEME: motif " + mprint(motif["plus"], motif["minus"]) + "; e-value exceeds cutoff")

        # motif(s) found :-)
        elif "Stopped because requested number of motifs (1) found" in reason:

            # find anchor genes' sequence_id
            training_set = e.findall("training_set/sequence") # all promoter sequences passed to MEME
            for i in xrange(len(training_set)):
                if "__ANCHOR" in training_set[i].attrib["name"]:
                    anchor_seq_id = training_set[i].attrib["id"] # e.g. id=sequence_1

            # only accept motifs which occur in the anchor genes promoter
            contributing_sites = e.findall("motifs/motif/contributing_sites/contributing_site") # sequences which contributed to the motif
            if anchor_seq_id in map(lambda site: site.attrib["sequence_id"], contributing_sites):
                # save motif score
                motif["score"] = e.find("motifs/motif").attrib["e_value"] # one motif, didn't ask MEME for more

                # save sequence sites which represent the motif
                motif["seqs"] = ["".join(map(lambda letter: letter.attrib["letter_id"], site.findall("site/letter_ref")))
                    for site in contributing_sites]

                # write sites to fasta file
                with open(os.path.join(meme_dir, mprint(motif["plus"], motif["minus"]), "binding_sites.fasta"), "w") as handle:
                    handle.write(">{}__{}\n".format(anchor, mprint(motif["plus"], motif["minus"])))
                    handle.write("\n".join(motif["seqs"]))

                logging.debug("MEME: motif {}; e-value = {}".format(
                    mprint(motif["plus"], motif["minus"]), motif["score"]))
            else:
                logging.debug("MEME: motif " + mprint(motif["plus"], motif["minus"]) + "; does not occur in anchor gene promoter")

        # unexpected reason, don't know why MEME stopped :-$
        else:
            logging.error("MEME stopped unexpectedly (reason: " + reason + ")")

    return filter(lambda m: m["score"] is not None, promoter_sets)


def filter_fimo_results(motifs, fimo_dir, promoters, anchor_promoter, options):
    """Analyse and filter FIMO results"""
    # TODO command-line option?
    if not hasattr(options, "max_percentage"):
        options.max_percentage = 14.0

    for motif in motifs:
        motif["hits"] = {}
        with open(os.path.join(fimo_dir, mprint(motif["plus"], motif["minus"]), "fimo.txt"), "r") as handle:
            table = csv.reader(handle, delimiter = "\t")
            for row in table:
                if not row[0].startswith("#"): # skip comment lines
                    seq_id = row[1]
                    if seq_id in motif["hits"]:
                        motif["hits"][seq_id] += 1
                    else:
                        motif["hits"][seq_id] = 1

        # write binding sites per promoter to file
        with open(os.path.join(fimo_dir, mprint(motif["plus"], motif["minus"]), "bs_per_promoter.csv"), "w") as handle:
            table = csv.writer(handle, delimiter = "\t", lineterminator = "\n")
            table.writerow(["#", "promoter", "binding sites"]) # table head
            for i in xrange(len(promoters)):
                promoter = get_promoter_id(promoters[i])
                if promoter in motif["hits"]:
                    table.writerow([i+1, promoter, motif["hits"][promoter]])
                else:
                    table.writerow([i+1, promoter, 0])

        percentage = float(len(motif["hits"])) / float(len(promoters)) * 100 # float!
        if percentage == 0.0:
            # too low
            logging.debug("FIMO: motif {}; occurs in {} promoters (no hits)".format(
                mprint(motif["plus"], motif["minus"]), len(motif["hits"])))
            motif["hits"] = None
        elif percentage > options.max_percentage:
            # too high
            logging.debug("FIMO: {}; occurs in {} promoters; {:.2f}% of all promoters (too many)".format(
                mprint(motif["plus"], motif["minus"]), len(motif["hits"]), percentage))
            motif["hits"] = None
        elif get_promoter_id(promoters[anchor_promoter]) not in motif["hits"]: # not in achor promoter
            # no site in anchor promoter
            logging.debug("FIMO: motif " + mprint(motif["plus"], motif["minus"]) + "; not hits in the promoter of the anchor gene")
            motif["hits"] = None
        else:
            # everything ok
            logging.debug("FIMO: motif {}; occurs in {} promoters; {:.2f}% of all promoters".format(
                mprint(motif["plus"], motif["minus"]), len(motif["hits"]), percentage))

    return filter(lambda m: m["hits"] is not None, motifs)


def get_islands(anchor_promoter, motifs, promoters, options):
    """Find islands of binding sites (previously found by FIMO) around anchor gene to define cluster borders"""
    max_gap_length = 2 # TODO options
    islands = []

    for motif in motifs:
        # create list with binding sites per promoter
        bs_per_promoter = [0] * len(promoters) # first: set number of binding sites to 0
        for i in xrange(len(promoters)):
            if get_promoter_id(promoters[i]) in motif["hits"]: # second: set actual number of binding sites, if any
                bs_per_promoter[i] = motif["hits"][get_promoter_id(promoters[i])]

        # upstream
        start = anchor_promoter # init upstream cluster border
        i = anchor_promoter # init position of anchor gene's promoter
        while i > 0:

            # promoter with binding site
            # … 1 …
            if bs_per_promoter[i-1] >= 1:
                start -= 1

            # no binding site, gap with length 1
            # … 1 0 1  …
            elif (i - 2 >= 0
                  and max_gap_length >= 1
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i-1] == 0
                  and bs_per_promoter[i-2] >= 1):
                start -= 2
                i -= 1

            # no binding site, gap with length 2
            # … 1 0 0 1  …
            elif (i - 3 >= 0
                  and max_gap_length >= 2
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i-1] == 0
                  and bs_per_promoter[i-2] == 0
                  and bs_per_promoter[i-3] >= 1):
                start -= 3
                i     -= 2

            # no binding site, gap with length 3
            # … 1 0 0 0 1  …
            elif (i - 4 >= 0
                  and max_gap_length >= 3
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i-1] == 0
                  and bs_per_promoter[i-2] == 0
                  and bs_per_promoter[i-3] == 0
                  and bs_per_promoter[i-4] >= 1):
                start -= 4
                i     -= 3

            # no binding site, gap with length 4
            # … 1 0 0 0 0 1  …
            elif (i - 5 >= 0
                  and max_gap_length >= 4
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i-1] == 0
                  and bs_per_promoter[i-2] == 0
                  and bs_per_promoter[i-3] == 0
                  and bs_per_promoter[i-4] == 0
                  and bs_per_promoter[i-5] >= 1):
                start -= 5
                i     -= 4

            # no binding site, gap with length 5
            # … 1 0 0 0 0 0 1  …
            elif (i - 6 >= 0
                  and max_gap_length >= 5
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i-1] == 0
                  and bs_per_promoter[i-2] == 0
                  and bs_per_promoter[i-3] == 0
                  and bs_per_promoter[i-4] == 0
                  and bs_per_promoter[i-5] == 0
                  and bs_per_promoter[i-6] >= 1):
                start -= 6
                i     -= 5

            # gap too long, stop upstream cluster extension
            else:
                break

            i -= 1

        # downstream
        i = anchor_promoter # reset position of anchor gene's promoter
        end = anchor_promoter # init downstream cluster border
        while i < len(bs_per_promoter) - 1:

            # promoter with binding site(s)
            # … 1 …
            if bs_per_promoter[i+1] > 0:
                end += 1

            # no binding site, gap with length 1
            # … 1 0 1  …
            elif (i + 2 < len(bs_per_promoter)
                  and max_gap_length >= 1
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i+1] == 0
                  and bs_per_promoter[i+2] >= 1):
                end += 2
                i += 1

            # no binding site, gap with length 2
            # … 1 0 0 1  …
            elif (i + 3 < len(bs_per_promoter)
                  and max_gap_length >= 2
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i+1] == 0
                  and bs_per_promoter[i+2] == 0
                  and bs_per_promoter[i+3] >= 1):
                end += 3
                i += 2

            # no binding site, gap with length 3
            # … 1 0 0 0 1  …
            elif (i + 4 < len(bs_per_promoter)
                  and max_gap_length >= 3
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i+1] == 0
                  and bs_per_promoter[i+2] == 0
                  and bs_per_promoter[i+3] == 0
                  and bs_per_promoter[i+4] >= 1):
                end += 4
                i += 3

            # no binding site, gap with length 4
            # … 1 0 0 0 0 1  …
            elif (i + 5 < len(bs_per_promoter)
                  and max_gap_length >= 4
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i+1] == 0
                  and bs_per_promoter[i+2] == 0
                  and bs_per_promoter[i+3] == 0
                  and bs_per_promoter[i+4] == 0
                  and bs_per_promoter[i+5] >= 1):
                end += 5
                i += 4

            # no binding site, gap with length 5
            # … 1 0 0 0 0 0 1  …
            elif (i + 6 < len(bs_per_promoter)
                  and max_gap_length >= 5
                  and bs_per_promoter[i] >= 1
                  and bs_per_promoter[i+1] == 0
                  and bs_per_promoter[i+2] == 0
                  and bs_per_promoter[i+3] == 0
                  and bs_per_promoter[i+4] == 0
                  and bs_per_promoter[i+5] == 0
                  and bs_per_promoter[i+6] >= 1):
                end += 6
                i += 5

            # gap too long, stop downstream cluster extension
            else:
                break

            i += 1

        logging.debug("Island {} -- {} (motif {})".format(
            get_promoter_id(promoters[start]), get_promoter_id(promoters[end]), mprint(motif["plus"], motif["minus"])))
        islands.append({"start": promoters[start], "end": promoters[end], "motif": motif})

    return islands


def sort_by_abundance(islands):
    """Sort upstream (start) and downstream (end) borders of islands by abundance"""
    # count border abundance
    # start/end are treated independently!
    starts = {}
    ends = {}
    for i in islands:
        if i["start"]["id"][0] in starts:
            starts[i["start"]["id"][0]]["abund"] += 1
        else:
            starts[i["start"]["id"][0]] = {"abund": 1}

        if i["end"]["id"][-1] in ends:
            ends[i["end"]["id"][-1]]["abund"] += 1
        else:
            ends[i["end"]["id"][-1]] = {"abund": 1}

        # keep track of motif score --> to sort by score if same abundance occurs more than once
        # AND
        # save motif "names" (plus, minus) --> additional info for showing the results later on
        if "mscore" not in starts[i["start"]["id"][0]] or float(i["motif"]["score"]) < float(starts[i["start"]["id"][0]]["mscore"]):
            starts[i["start"]["id"][0]]["mscore"] = i["motif"]["score"]
            starts[i["start"]["id"][0]]["plus"] = i["motif"]["plus"]
            starts[i["start"]["id"][0]]["minus"] = i["motif"]["minus"]
        if "mscore" not in ends[i["end"]["id"][-1]] or float(i["motif"]["score"]) < float(ends[i["end"]["id"][-1]]["mscore"]):
            ends[i["end"]["id"][-1]]["mscore"] = i["motif"]["score"]
            ends[i["end"]["id"][-1]]["plus"] = i["motif"]["plus"]
            ends[i["end"]["id"][-1]]["minus"] = i["motif"]["minus"]

    # compute sum of start and end abundance, remove duplicates, sort descending
    abundances_sum_sorted = sorted(set([s["abund"] + e["abund"] for s in starts.values() for e in ends.values()]), reverse=True)
    # compute sum of start and end motif score, remove duplicates, sort ascending
    scores_sum_sorted = sorted(set([float(s["mscore"]) + float(e["mscore"]) for s in starts.values() for e in ends.values()]))
    # sort by value (=abundance) of start, descending
    starts_sorted = sorted(starts, key=starts.get, reverse=True)
    # sort by value (=abundance) of end, descending
    ends_sorted = sorted(ends, key=ends.get, reverse=True)

    clusters = []
    for abundance in abundances_sum_sorted:
        # list from highest (best) to lowest (worst) abundance
        for score in scores_sum_sorted:
            # list from lowest (best) to highest (worst) motif score/e-value
            for start in starts_sorted:
                for end in ends_sorted:
                    if (starts[start]["abund"] + ends[end]["abund"] == abundance
                            and float(starts[start]["mscore"]) + float(ends[end]["mscore"]) == score):
                        clusters.append({
                            "start": {
                                "gene": start,
                                "abundance": starts[start]["abund"],
                                "score": starts[start]["mscore"],
                                "plus": starts[start]["plus"],
                                "minus": starts[start]["minus"],
                            },
                            "end": {
                                "gene": end,
                                "abundance": ends[end]["abund"],
                                "score": ends[end]["mscore"],
                                "plus": ends[end]["plus"],
                                "minus": ends[end]["minus"],
                            },
                        })

                        logging.debug("Upstream border:   gene {}; abundance {}; motif {}; score {}".format(
                            start,
                            starts[start]["abund"],
                            mprint(starts[start]["plus"], starts[start]["minus"]),
                            starts[start]["mscore"],
                        ))
                        logging.debug("Downstream border: gene {}; abundance {}; motif {}; score {}".format(
                            end,
                            ends[end]["abund"],
                            mprint(ends[end]["plus"], ends[end]["minus"]),
                            ends[end]["mscore"],
                        ))
                        logging.debug("Total abundance {}, total score {:.1e}".format(
                            abundance, score))

    return clusters


def check_cluster_predictions(cluster_predictions, seq_record, promoters, ignored_genes):
    """Get some more infos about each cluster prediction and check if it seems to be sane"""
    checked_predictions = []
    for cp in xrange(len(cluster_predictions)):
        prediction = cluster_predictions[cp]
        sane = True

        # find indices of first and last GENE of the cluster prediction in all genes
        all_genes = map(lambda g: utils.get_gene_id(g), utils.get_all_features_of_type(seq_record, "gene"))

        start_index_genes = None
        end_index_genes = None
        for i in xrange(len(all_genes)):
            if not start_index_genes and prediction["start"]["gene"] == all_genes[i]:
                start_index_genes = i
            if not end_index_genes and prediction["end"]["gene"] == all_genes[i]:
                end_index_genes = i
            if start_index_genes and end_index_genes:
                break

        # find indices of first and last PROMOTER of the cluster prediction in all promoters
        start_index_promoters = None
        end_index_promoters = None
        for i in xrange(len(promoters)):
            if not start_index_promoters and prediction["start"]["gene"] in promoters[i]["id"]:
                start_index_promoters = i
            if not end_index_promoters and prediction["end"]["gene"] in promoters[i]["id"]:
                end_index_promoters = i
            if start_index_promoters and end_index_promoters:
                break

        prediction["start"]["promoter"] = get_promoter_id(promoters[start_index_promoters])
        prediction["end"]["promoter"] = get_promoter_id(promoters[end_index_promoters])
        prediction["genes"] = end_index_genes - start_index_genes + 1
        prediction["promoters"] = end_index_promoters - start_index_promoters + 1

        if cp == 0:
            logging.debug("Best prediction (most abundant): {!r} -- {!r}".format(
                prediction["start"]["gene"], prediction["end"]["gene"]))
        else:
            logging.debug("Alternative prediction ({}): {!r} -- {!r}".format(
                cp, prediction["start"]["gene"], prediction["end"]["gene"]))

        # warn if cluster prediction right at or next to record (~ contig) border
        if start_index_genes < 10:
            logging.debug(
                "Upstream cluster border located at or next to sequence record border, prediction could have been truncated by record border")
            sane = False
        if end_index_genes > len(all_genes) - 10:
            logging.debug(
                "Downstream cluster border located at or next to sequence record border, prediction could have been truncated by record border")
            sane = False

        # warn if cluster prediction too short (includes less than 3 genes)
        if prediction["genes"] < 3:
            logging.debug("Cluster is very short (less than 3 genes). Prediction may be questionable.")
            sane = False

        # warn if ignored gene (overlapping with anthor gene, see ignore_overlapping()) would have been part of the cluster
        for ignored_gene in map(lambda g: utils.get_gene_id(g), ignored_genes):
            if ignored_gene in all_genes[start_index_genes : end_index_genes + 1]:
                logging.debug("Gene {!r} is part of the predicted cluster, but it is overlapping with another gene and was ignored".format(
                    ignored_gene))
                logging.debug("Gene {!r} could have effected the cluster prediction".format(
                    ignored_gene))
                # sane = False # uncomment if you want alternatives for predictions with ignored genes, too
                break

        checked_predictions.append(prediction)

        if sane:
            break

    return checked_predictions


def cleanup_outdir(anchor_genes, cluster_predictions, options):
    """Delete unnecessary files to free disk space"""
    all_motifs = set()
    for motif in _plus_minus:
        all_motifs.add(mprint(motif["plus"], motif["minus"]))

    for anchor in anchor_genes:
        if anchor in cluster_predictions:
            used_motifs = set()
            for cluster in cluster_predictions[anchor]:
                used_motifs.add(mprint(cluster["start"]["plus"], cluster["start"]["minus"]))
                used_motifs.add(mprint(cluster["end"]["plus"], cluster["end"]["minus"]))
            unused_motifs = all_motifs.difference(used_motifs)
            # only remove directories from "unused" motifs (no cluster prediction)
            for directory in unused_motifs:
                shutil.rmtree(os.path.join(options.outputfoldername, "meme", anchor, directory), ignore_errors=True)
                shutil.rmtree(os.path.join(options.outputfoldername, "fimo", anchor, directory), ignore_errors=True)
        else:
            # all motifs are "unused" (not a single prediction for this anchor gene)
            # --> remove anchor genes directory, including all motif subdirectories
            shutil.rmtree(os.path.join(options.outputfoldername, "meme", anchor), ignore_errors=True)
            shutil.rmtree(os.path.join(options.outputfoldername, "fimo", anchor), ignore_errors=True)


### storage methods ###
def store_promoters(promoters, seq_record):
    """Store information about promoter sequences to a SeqRecord"""
    for promoter in promoters:
        new_feature = SeqFeature.SeqFeature(
            # 1-based GenBank vs. 0-based Python format --> -1 is important!
            # see http://biopython.org/DIST/docs/api/Bio.SeqFeature.FeatureLocation-class.html
            FeatureLocation(promoter["start"] - 1, promoter["end"]), type = "promoter")
        new_feature.qualifiers = {
            "locus_tag" : promoter["id"], # already a list with one or two elements
            "seq"       : [str(promoter["seq"])], # TODO save string or Seq object?
        }

        if len(promoter["id"]) > 1:
            new_feature.qualifiers["note"] = ["bidirectional promoter"]

        seq_record.features.append(new_feature)


def store_clusters(anchor, clusters, seq_record):
    """Store the borders of predicted clusters to a SeqRecord"""
    cluster_features = utils.get_cluster_features(seq_record)
    for i in xrange(len(clusters)):
        cluster = clusters[i]
        # cluster borders returned by hmmdetect are based on CDS features
        # see find_clusters() in hmmdetect's __init__.py
        # in contrast, cluster borders returned by cassis are based on gene features
        # --> hmmdetect derived clusters have exact loctions, like the CDSs have
        # --> cassis derived clusters may have fuzzy locations, like the genes have
        #
        # utils.get_all_features_of_type_with_query() returns a list
        # there should be no second gene with the same locus tag
        # --> always take the first [0] element of the return value
        left = utils.get_all_features_of_type_with_gene_id(
            seq_record, "gene", cluster["start"]["gene"])[0]
        right = utils.get_all_features_of_type_with_gene_id(
            seq_record, "gene", cluster["end"]["gene"])[0]

        new_feature = SeqFeature.SeqFeature(
            FeatureLocation(left.location.start, right.location.end), type = "cluster_border")
        new_feature.qualifiers = {
            "tool"             : ["cassis"],
            "anchor"           : [anchor],
            "abundance"        : [cluster["start"]["abundance"] + cluster["end"]["abundance"]],
            "motif_score"      : ["{:.1e}".format(float(cluster["start"]["score"]) + float(cluster["end"]["score"]))],
            "gene_left"        : [cluster["start"]["gene"]],
            "promoter_left"    : [cluster["start"]["promoter"]],
            "abundance_left"   : [cluster["start"]["abundance"]],
            "motif_left"       : [mprint(cluster["start"]["plus"], cluster["start"]["minus"])],
            "motif_score_left" : ["{:.1e}".format(float(cluster["start"]["score"]))],
            "gene_right"       : [cluster["end"]["gene"]],
            "promoter_right"   : [cluster["end"]["promoter"]],
            "abundance_right"  : [cluster["end"]["abundance"]],
            "motif_right"      : [mprint(cluster["end"]["plus"], cluster["end"]["minus"])],
            "motif_score_right": ["{:.1e}".format(float(cluster["end"]["score"]))],
            "genes"            : [cluster["genes"]],
            "promoters"        : [cluster["promoters"]],
        }

        if i == 0:
            new_feature.qualifiers["note"] = ["best prediction (most abundant) for anchor gene {}".format(anchor)]
        else:
            new_feature.qualifiers["note"] = ["alternative prediction ({}) for anchor gene {}".format(i, anchor)]

        seq_record.features.append(new_feature)

        # only extend existing clusters for the best prediction
        if i != 0:
            continue

        for existing_cluster in cluster_features:
            if not utils.features_overlap(new_feature, existing_cluster):
                continue

            # extend existing cluster start location
            if existing_cluster.location.start > new_feature.location.start:
                existing_cluster.location = FeatureLocation(int(new_feature.location.start),
                                                            existing_cluster.location.end)

            # extend cluster end location
            if existing_cluster.location.end < new_feature.location.end:
                existing_cluster.location = FeatureLocation(existing_cluster.location.start,
                                                            int(new_feature.location.end))
