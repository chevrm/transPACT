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
import codecs
from collections import OrderedDict
import logging
import os
import sys
import shutil
from os import path
import subprocess
import string
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
from argparse import Namespace
import warnings
# Don't display the SearchIO experimental warning, we know this.
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO
from Bio import SeqIO
from Bio.Alphabet import generic_protein, generic_dna
from Bio.Seq import UnknownSeq
from antismash.config import get_config
from helperlibs.wrappers.io import TemporaryDirectory
import signal
from multiprocessing import Pool

from zipfile import ZipFile, ZIP_DEFLATED, LargeZipFile
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import re
import argparse

from numpy import array_split, array
import urllib2
from urllib2 import URLError
import httplib
import time


class Storage(dict):
    """Simple storage class"""
    def __init__(self, indict=None):
        if indict is None:
            indict = {}
        dict.__init__(self, indict)
        self.__initialized = True

    def __getattr__(self, attr):
        try:
            return self.__getitem__(attr)
        except KeyError:
            raise AttributeError(attr)

    def __setattr__(self, attr, value):
        if '_Storage__initialized' not in self.__dict__:
            return dict.__setattr__(self, attr, value)
        elif attr in self:
            dict.__setattr__(self, attr, value)
        else:
            self.__setitem__(attr, value)


def getArgParser():
    class AntiSmashParser(argparse.ArgumentParser):
        """Custom argument parser for antiSMASH
        """
        _showAll = False
        _displayGroup = {}

        def __init__(self, *args):
            """Initialisation method for the parser class"""
            kwargs = {}
            kwargs["add_help"] = False
            super(AntiSmashParser, self).__init__(*args, **kwargs)

        def add_argument_group(self, *args, **kwargs):
            if not args[0] in self._displayGroup:
                self._displayGroup[args[0]] = []
            if "basic" in kwargs:
                if kwargs["basic"]:
                    self._displayGroup[args[0]].extend(["basic"])
                del kwargs["basic"]
            if "param" in kwargs:
                self._displayGroup[args[0]].extend(kwargs["param"])
                del kwargs["param"]
            group = super(AntiSmashParser, self).add_argument_group(*args, **kwargs)
            return group

        def print_help(self, file=None, showAll=False):
            self._showAll = showAll
            super(AntiSmashParser, self).print_help(file)

        def format_help(self):
            """Custom help format"""
            help_text = """
########### antiSMASH ver. {version} #############

{usage}

{args}
--------
Options
--------
{opts}
""".format(version=get_version(), usage=self.format_usage(), args=self._get_args_text(), opts=self._get_opts_text())
            return help_text

        def format_usage(self):
            if self._showAll:
                formatter = self._get_formatter()
                formatter.add_usage(self.usage, self._actions, self._mutually_exclusive_groups)
                return formatter.format_help()
            return "usage: {prog} [-h] [options ..] [sequence [sequence ..]]".format(prog=self.prog) + "\n"

        def _get_args_text(self):
            # fetch arg lists using formatter
            formatter = self._get_formatter()
            for action_group in self._action_groups:
                if action_group.title == "positional arguments":
                    formatter.start_section("arguments")
                    formatter.add_arguments(action_group._group_actions)
                    formatter.end_section()
                    break
            return formatter.format_help()

        def _get_opts_text(self):
            # fetch opt lists using formatter
            formatter = self._get_formatter()
            for action_group in self._action_groups:
                if action_group.title in ["optional arguments"]:
                    formatter.add_arguments(action_group._group_actions)
            for action_group in self._action_groups:
                if action_group.title not in ["optional arguments", "positional arguments"]:
                    show_opt = self._showAll
                    if (not show_opt):
                        if "basic" in self._displayGroup[action_group.title]:
                            show_opt = True
                        elif len(list(set(sys.argv) & set(self._displayGroup[action_group.title]))) > 0:
                            show_opt = True
                    if show_opt:
                        formatter.start_section(action_group.title)
                        if action_group.description is None:
                            action_group.description = ''
                        formatter.add_text(action_group.description)
                        formatter.add_arguments(action_group._group_actions)
                        formatter.end_section()
            return formatter.format_help()

    return AntiSmashParser()


def get_all_features_of_type(seq_record, types):
    "Return all features of the specified types for a seq_record"
    if isinstance(types, str):
        # force into a tuple
        types = (types, )
    features = []
    for f in seq_record.features:
        if f.type in types:
            features.append(f)
    return features

def get_all_features_of_type_with_query(seq_record, feature_type, query_tag, query_value):
    """Return all features of type 'type' which contain a 'query_tag' with value 'query_value'
    Note: query has to be exact!"""

    features_with_type = get_all_features_of_type(seq_record, feature_type)
    features = []
    for feature_to_test in features_with_type:

        if (query_tag in feature_to_test.qualifiers) and (query_value in feature_to_test.qualifiers[query_tag]):
            features.append(feature_to_test)
    return features


def get_all_features_of_type_with_gene_id(seq_record, feature_type, gene_id):
    """Return all features of the given type with a matching gene_id"""
    features_with_type = get_all_features_of_type(seq_record, feature_type)
    features = []
    for feature_to_test in features_with_type:
        if get_gene_id(feature_to_test) == gene_id:
            features.append(feature_to_test)
    return features


def get_cds_features(seq_record):
    "Return all CDS features for a seq_record"
    return get_all_features_of_type(seq_record, "CDS")


def get_cluster_border_features(seq_record):
    """Return all cluster_border features for a seq_record"""
    return get_all_features_of_type(seq_record, "cluster_border")


def get_withincluster_cds_features(seq_record):
    features = get_cds_features(seq_record)
    clusters = get_cluster_features(seq_record)
    withinclusterfeatures = []
    for feature in features:
        for cluster in clusters:
            if not (cluster.location.start <= feature.location.start <= cluster.location.end or \
               cluster.location.start <= feature.location.end <= cluster.location.end):
                continue
            if feature not in withinclusterfeatures:
                withinclusterfeatures.append(feature)
    return withinclusterfeatures

def get_cluster_cds_features(cluster, seq_record):
    clustercdsfeatures = []
    for feature in seq_record.features:
        if feature.type != 'CDS':
            continue
        if cluster.location.start <= feature.location.start <= cluster.location.end or \
           cluster.location.start <= feature.location.end <= cluster.location.end:
            clustercdsfeatures.append(feature)
    return clustercdsfeatures


def get_cluster_cluster_border_features(cluster, seq_record):
    """Get cluster_border features that overlap with the given cluster"""
    borders = get_cluster_border_features(seq_record)
    overlapping_borders = []
    for border in borders:
        if features_overlap(cluster, border):
            overlapping_borders.append(border)

    return overlapping_borders


def get_cluster_aSDomain_features(cluster, seq_record):
    aSDomainfeatures = []
    for feature in seq_record.features:
        if feature.type != 'aSDomain':
            continue
        if cluster.location.start <= feature.location.start <= cluster.location.end or \
           cluster.location.start <= feature.location.end <= cluster.location.end:
            aSDomainfeatures.append(feature)
    return aSDomainfeatures

def features_overlap(a, b):
    "Check if two features have overlapping locations"
    astart = min([a.location.start, a.location.end])
    aend = max([a.location.start, a.location.end])
    bstart = min([b.location.start, b.location.end])
    bend = max([b.location.start, b.location.end])
    return (astart >= bstart and astart <= bend) or \
           (aend >= bstart and aend <= bend) or \
           (bstart >= astart and bstart <= aend) or \
           (bend >= astart and bend <= aend)

def get_secmet_cds_features(seq_record):
    features = get_cds_features(seq_record)
    secmet_features = []
    for feature in features:
        if 'sec_met' in feature.qualifiers:
            for annotation in feature.qualifiers['sec_met']:
                if annotation.startswith('Type: ') and not annotation.endswith('none'):
                    secmet_features.append(feature)
                    break
    return secmet_features

def get_pksnrps_cds_features(seq_record):
    features = get_cds_features(seq_record)
    pksnrpscoregenes = []
    for feature in features:
        if 'sec_met' in feature.qualifiers:
            for annotation in feature.qualifiers['sec_met']:
                if annotation.startswith('NRPS/PKS Domain:'):
                    pksnrpscoregenes.append(feature)
                    break
    return pksnrpscoregenes

def get_nrpspks_domain_dict(seq_record):
    domaindict = {}
    features = get_cds_features(seq_record)
    for feature in features:
        domainlist = []
        if 'sec_met' in feature.qualifiers:
            domains = [qualifier for qualifier in feature.qualifiers['sec_met'] if "NRPS/PKS Domain: " in qualifier]
            for domain in domains:
                hit_id =  domain.partition("NRPS/PKS Domain: ")[2].partition(" (")[0]
                domstart = domain.partition("NRPS/PKS Domain: ")[2].partition(" (")[2].partition("-")[0]
                domend = domain.partition("NRPS/PKS Domain: ")[2].partition("). ")[0].rpartition("-")[2]
                evalue = domain.partition("E-value: ")[2].partition(". Score:")[0]
                bitscore = domain.partition("Score: ")[2].partition(";")[0]
                domainlist.append([hit_id, int(domstart), int(domend), evalue, float(bitscore)])
            if len(domainlist) > 0:
                domaindict[get_gene_id(feature)] = domainlist
    return domaindict

def get_nrpspks_substr_spec_preds(seq_record):
    substr_spec_preds = Storage()
    substr_spec_preds.consensuspreds = {}
    substr_spec_preds.nrps_svm_preds = {}
    substr_spec_preds.nrps_code_preds = {}
    substr_spec_preds.minowa_nrps_preds = {}
    substr_spec_preds.pks_code_preds = {}
    substr_spec_preds.minowa_pks_preds = {}
    substr_spec_preds.minowa_cal_preds = {}
    substr_spec_preds.kr_activity_preds = {}
    substr_spec_preds.kr_stereo_preds = {}
    features = get_cds_features(seq_record)
    for feature in features:
        nrat, nra, nrcal, nrkr = 0, 0, 0, 0
        if 'sec_met' in feature.qualifiers:
            domains = [qualifier for qualifier in feature.qualifiers['sec_met'] if "NRPS/PKS Domain: " in qualifier]
            for domain in domains:
                if "AMP-binding" in domain or "A-OX" in domain:
                    nra += 1
                    domainname = get_gene_id(feature) + "_A" + str(nra)
                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
                    nrps_svm_pred = predictionstext.partition(" (NRPSPredictor2 SVM)")[0]
                    nrps_code_pred = predictionstext.partition(" (NRPSPredictor2 SVM), ")[2].partition(" (Stachelhaus code)")[0]
                    minowa_nrps_pred = predictionstext.partition("(Stachelhaus code), ")[2].partition(" (Minowa)")[0]
                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
                    substr_spec_preds.nrps_svm_preds[domainname] = nrps_svm_pred
                    substr_spec_preds.nrps_code_preds[domainname] = nrps_code_pred
                    substr_spec_preds.minowa_nrps_preds[domainname] = minowa_nrps_pred
                    substr_spec_preds.consensuspreds[domainname] = consensuspred
                elif "PKS_AT" in domain:
                    nrat += 1
                    domainname = get_gene_id(feature) + "_AT" + str(nrat)
                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
                    pks_code_pred = predictionstext.partition(" (PKS signature)")[0]
                    minowa_pks_pred = predictionstext.partition("(PKS signature), ")[2].partition(" (Minowa)")[0]
                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
                    substr_spec_preds.pks_code_preds[domainname] = pks_code_pred
                    substr_spec_preds.minowa_pks_preds[domainname] = minowa_pks_pred
                    substr_spec_preds.consensuspreds[domainname] = consensuspred
                elif "CAL_domain" in domain:
                    nrcal += 1
                    domainname = get_gene_id(feature) + "_CAL" + str(nrcal)
                    predictionstext = domain.partition("Substrate specificity predictions:")[2]
                    minowa_cal_pred = predictionstext.partition(" (Minowa)")[0]
                    consensuspred = predictionstext.partition("(Minowa), ")[2].partition(" (consensus)")[0]
                    substr_spec_preds.minowa_cal_preds[domainname] = minowa_cal_pred
                    substr_spec_preds.consensuspreds[domainname] = consensuspred
                elif "PKS_KR" in domain:
                    nrkr += 1
                    domainname = get_gene_id(feature) + "_KR" + str(nrkr)
                    activityprediction = domain.partition("Predicted KR activity: ")[2].partition(";")[0]
                    stereoprediction = domain.partition("Predicted KR stereochemistry: ")[2].partition(";")[0]
                    substr_spec_preds.kr_activity_preds[domainname] = activityprediction
                    substr_spec_preds.kr_stereo_preds[domainname] = stereoprediction
    return substr_spec_preds

def get_smcog_annotations(seq_record):
    smcogdict = {}
    smcogdescriptions = {}
    features = get_cds_features(seq_record)
    for feature in features:
        if 'note' in feature.qualifiers:
            notes = feature.qualifiers['note']
            for note in notes:
                if "smCOG: " in note:
                    smcogid = note.partition("smCOG: ")[2].partition(":")[0]
                    smcog_descr = note.partition("smCOG: ")[2].partition(":")[2].partition("(Score:")[0]
                    smcogdict[get_gene_id(feature)] = smcogid
                    smcogdescriptions[smcogid] = smcog_descr
    return smcogdict, smcogdescriptions

def get_pfam_features(seq_record):
    "Return all CDS_motif features containing a 'PFAM-Id: ' note for a seq_record"
    pfam_features = []
    for feature in get_all_features_of_type(seq_record, "PFAM_domain"):
        if 'db_xref' not in feature.qualifiers:
            continue
        for xref in feature.qualifiers['db_xref']:
            if xref.startswith("PFAM: "):
                pfam_features.append(feature)
                break
    return pfam_features

def get_cluster_features(seq_record):
    "Return all cluster features for a seq_record"
    return get_all_features_of_type(seq_record, "cluster")

def get_sorted_cluster_features(seq_record):
    "Return all cluster features for a seq_record"
    clusters = get_all_features_of_type(seq_record, "cluster")
    if len(clusters) == 0:
        return []
    numberdict = {}
    for cluster in clusters:
        numberdict[get_cluster_number(cluster)] = cluster
    return [numberdict[clusternr] for clusternr in numberdict.keys()]

def get_structure_pred(cluster):
    "Return all structure prediction for a cluster feature"
    if 'note' in cluster.qualifiers:
        for note in cluster.qualifiers['note']:
            if "Monomers prediction: " in note:
                return note.partition("Monomers prediction: ")[2]
    if get_cluster_type(cluster) == 'ectoine':
        return 'ectoine'
    return "N/A"

def get_cluster_number(cluster):
    "Get the integer representation of the Cluster number qualifier"
    if 'note' in cluster.qualifiers and "Cluster number: " in cluster.qualifiers['note'][0]:
        clusternr = int(cluster.qualifiers['note'][0].partition("Cluster number: ")[2])
        return clusternr
    else:
        return 0

def get_cluster_type(cluster):
    "Get product type of a gene cluster"
    return cluster.qualifiers['product'][0]

def get_cluster_by_nr(seq_record, queryclusternr):
    "Return all cluster features for a seq_record of a certain type"
    clusters = get_all_features_of_type(seq_record, "cluster")
    for cluster in clusters:
        if "Cluster number: " in cluster.qualifiers['note'][0]:
            clusternr = int(cluster.qualifiers['note'][0].partition("Cluster number: ")[2])
            if clusternr == queryclusternr:
                return cluster

def get_cluster_features_of_type(seq_record, clustertype):
    "Return all cluster features for a seq_record of a certain type"
    clusters = get_all_features_of_type(seq_record, "cluster")
    return [cluster for cluster in clusters if clustertype in cluster.qualifiers['product'][0]]

def locate_executable(name):
    "Find an executable in the path and return the full path"
    # In windows, executables tend to end on .exe
    if sys.platform == 'win32':
        name += ".exe"
    file_path, _ = os.path.split(name)
    if file_path != "":
        if path.isfile(name) and os.access(name, os.X_OK):
            logging.debug("Found executable %r", name)
            return name
    for p in os.environ["PATH"].split(os.pathsep):
        full_name = path.join(p, name)
        if path.isfile(full_name) and os.access(full_name, os.X_OK):
            logging.debug("Found executable %r", full_name)
            return full_name

    return None

def locate_file(name):
    "Find a file and return the full path"
    file_path, _ = os.path.split(name)
    if file_path != "":
        if path.isfile(name) and os.access(name, os.R_OK):
            logging.debug("Found file %r", name)
            return name
    return None

# Ignore the pylint warning about input being redifined, as we're just
# following the subprocess names here.
# pylint: disable=redefined-builtin
def execute(commands, input=None):
    "Execute commands in a system-independent manner"

    if input is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    try:
        proc = subprocess.Popen(commands, stdin=stdin_redir,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        out, err = proc.communicate(input=input)
        retcode = proc.returncode
        return out, err, retcode
    except OSError as e:
        logging.debug("%r %r returned %r", commands, input[:40] if input is not None else None, e)
        raise
# pylint: enable=redefined-builtin


def run_p(options, args):
    """Execute commands in a parallel manner"""
    os.setpgid(0, 0)
    p = Pool(options.cpus)
    jobs = p.map_async(child_process, args)

    try:
        errors = jobs.get()
    except KeyboardInterrupt:
        logging.error("Interrupted by user")
        os.killpg(os.getpid(), signal.SIGTERM)
        return 128

    p.close()

    return errors


def child_process(args):
    """Called by multiprocessing's map or map_async method"""
    try:
        logging.debug("Calling {!r}".format(" ".join(args)))
        p = subprocess.Popen(args) # actually only a chunk of args
        p.communicate()
    except KeyboardInterrupt:
        #  Need to raise some runtime error that is not KeyboardInterrupt, because Python has a bug there
        raise RuntimeError("Killed by keyboard interrupt")

    return p.returncode

def run_fimo_simple(query_motif_file, target_sequence):
    "Run FIMO"
    command = ["fimo", "--text", "--verbosity", "1", query_motif_file, target_sequence]
    try:
        out, err, retcode = execute(command)
    except OSError:
        return []
    if retcode != 0:
        logging.debug('FIMO returned %d: %r while searching %r', retcode,
                        err, query_motif_file)
        return []
    return out


def run_hmmsearch(query_hmmfile, target_sequence, use_tempfile=False):
    "Run hmmsearch"
    config = get_config()
    command = ["hmmsearch", "--cpu", str(config.cpus),
               "-o", os.devnull, # throw away the verbose output
               "--domtblout", "result.domtab",
               query_hmmfile]

    # Allow to disable multithreading for HMMer3 calls in the command line
    if 'hmmer3' in config and 'multithreading' in config.hmmer3 and \
            not config.hmmer3.multithreading:
        command = command[0:1] + command[3:]

    with TemporaryDirectory(change=True):
        try:
            if use_tempfile:
                with open("input.fa", 'w') as handle:
                    handle.write(target_sequence)
                command.append("input.fa")
                _, err, retcode = execute(command)
            else:
                command.append('-')
                _, err, retcode = execute(command, input=target_sequence)
        except OSError:
            return []
        if retcode != 0:
            logging.debug('hmmsearch returned %d: %r while searching %r', retcode,
                          err, query_hmmfile)
            return []
        results = list(SearchIO.parse("result.domtab", 'hmmsearch3-domtab'))
        return results


def run_hmmscan(target_hmmfile, query_sequence, opts=None, results_file=None):
    "Run hmmscan"
    config = get_config()
    command = ["hmmscan", "--cpu", str(config.cpus), "--nobias"]

    # Allow to disable multithreading for HMMer3 calls in the command line
    if 'hmmer3' in config and 'multithreading' in config.hmmer3 and \
            not config.hmmer3.multithreading:
        command = command[0:1] + command[3:]

    if opts is not None:
        command.extend(opts)
    command.extend([target_hmmfile, '-'])
    try:
        out, err, retcode = execute(command, input=query_sequence)
    except OSError:
        return []
    if retcode != 0:
        logging.debug('hmmscan returned %d: %r while scanning %r' , retcode,
                        err, query_sequence[:100])
        return []

    if results_file is not None:
        with open(results_file, 'w') as fh:
            fh.write(out)

    res_stream = StringIO(out)
    results = list(SearchIO.parse(res_stream, 'hmmer3-text'))
    return results


def run_hmmpfam2(query_hmmfile, target_sequence):
    "Run hmmpfam2"
    config = get_config()
    command = ["hmmpfam2", "--cpu", str(config.cpus),
               query_hmmfile, '-']

    # Allow to disable multithreading for HMMer2 calls in the command line
    if 'hmmer2' in config and 'multithreading' in config.hmmer2 and \
            not config.hmmer2.multithreading:
        command = command[0:1] + command[3:]

    try:
        out, err, retcode = execute(command, input=target_sequence)
    except OSError:
        return []
    if retcode != 0:
        logging.debug('hmmpfam2 returned %d: %r while searching %r', retcode,
                        err, query_hmmfile)
        return []
    res_stream = StringIO(out)
    results = list(SearchIO.parse(res_stream, 'hmmer2-text'))
    return results


def run_hmmpress(hmmfile):
    "Run hmmpress"
    command = ['hmmpress', hmmfile]
    try:
        out, err, retcode = execute(command)
    except OSError as e:
        retcode = 1
        err = str(e)
        out = None
    return out, err,  retcode


def run_meme(meme_dir, options):
    """Set paths, check existing files and run MEME in parallel on each promoter set"""
    args = []
    for plus_minus in os.listdir(meme_dir):
        input_file = os.path.join(meme_dir, plus_minus, "promoters.fasta")
        output_file = os.path.join(meme_dir, plus_minus, "meme.xml")

        # input file present and size not zero --> should be fine to use
        if os.path.isfile(input_file) and os.path.getsize(input_file) > 0:
            # output file already present and size not zero --> do not run MEME again on this one
            if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
                pass
            else:
                args.append([
                    "meme", input_file,
                    "-oc", os.path.dirname(input_file),
                    "-dna",
                    "-nostatus",
                    "-mod", "anr",
                    "-nmotifs", "1",
                    "-minw", "6",
                    "-maxw", "12",
                    "-revcomp",
                    "-evt", "1.0e+005",
                ])

    errors = run_p(options, args)
    return sum(errors)


def run_fimo(meme_dir, fimo_dir, seq_record, options):
    """Set paths, check existing files and run FIMO in parallel on each predicted motif"""
    if not os.path.exists(fimo_dir):
        os.makedirs(fimo_dir)

    args = []
    for plus_minus in os.listdir(meme_dir):
        motif_file = os.path.join(meme_dir, plus_minus, "meme.html")
        sites_file = os.path.join(meme_dir, plus_minus, "binding_sites.fasta")
        output_file = os.path.join(fimo_dir, plus_minus, "fimo.txt")
        output_dir = os.path.join(fimo_dir, plus_minus)

        # input file and binding sites file present and size not zero --> should be fine to use
        if (os.path.isfile(motif_file)
                and os.path.getsize(motif_file) > 0
                and os.path.isfile(sites_file)
                and os.path.getsize(sites_file) > 0):
            # output file already present and size not zero --> do not run FIMO again on this one
            if os.path.isfile(output_file) and os.path.getsize(output_file) > 0:
                pass
            else:
                args.append([
                    "fimo",
                    "-verbosity", "1",
                    "-motif", "1",
                    "-thresh", "0.00006",
                    "-oc", output_dir,
                    motif_file,
                    os.path.join(options.outputfoldername, seq_record.name + "_promoter_sequences.fasta"),
                ])

    errors = run_p(options, args)
    return sum(errors)


def hmmlengths(hmmfile):
    hmmlengthsdict = {}
    openedhmmfile = open(hmmfile,"r")
    filetext = openedhmmfile.read()
    filetext = filetext.replace("\r","\n")
    hmms = filetext.split("//")[:-1]
    for i in hmms:
        namepart = i.split("NAME  ")[1]
        name = namepart.split("\n")[0]
        lengthpart = i.split("LENG  ")[1]
        length = lengthpart.split("\n")[0]
        hmmlengthsdict[name] = int(length)
    return hmmlengthsdict

def cmp_feature_location(a, b):
    "Compare two features by their start/end locations"
    ret = cmp(a.location.start, b.location.start)
    if ret != 0:
        return ret
    return cmp(a.location.end, b.location.end)

def sort_features(seq_record):
    "Sort features in a seq_record by their position"
    #Check if all features have a proper location assigned
    for feature in seq_record.features:
        if feature.location is None:
            if feature.id != "<unknown id>":
                logging.error("Feature '%s' has no proper location assigned", feature.id)
            elif "locus_tag" in feature.qualifiers:
                logging.error("Feature '%s' has no proper location assigned", feature.qualifiers["locus_tag"][0])
            else:
                logging.error("File contains feature without proper location assignment")
            sys.exit(0) #FIXME: is sys.exit(0) really what we want to do here?
    #Sort features by location
    seq_record.features.sort(cmp=cmp_feature_location)

def fix_locus_tags(seq_record, config):
    "Fix CDS feature that don't have a locus_tag, gene name or protein id"
    if 'next_locus_tag' not in config:
        config.next_locus_tag = 1

    cds_list = get_cds_features(seq_record)
    for feature in cds_list:
        if get_gene_id(feature) == "no_tag_found":
            feature.qualifiers['locus_tag'] = ['AUTOORF_%05d' % config.next_locus_tag]
            config.next_locus_tag += 1
        #Fix locus tags, gene names or protein IDs if they contain illegal chars
        illegal_chars  = '''!"#$%&()*+,:; \r\n\t=>?@[]^`'{|}/ '''
        if 'locus_tag' in feature.qualifiers:
            for char in feature.qualifiers['locus_tag'][0]:
                if char in illegal_chars:
                    feature.qualifiers['locus_tag'][0] = feature.qualifiers['locus_tag'][0].replace(char, "_")
        if 'gene' in feature.qualifiers:
            for char in feature.qualifiers['gene'][0]:
                if char in illegal_chars:
                    feature.qualifiers['gene'][0] = feature.qualifiers['gene'][0].replace(char, "_")
        if 'protein_id' in feature.qualifiers:
            for char in feature.qualifiers['protein_id'][0]:
                if char in illegal_chars:
                    feature.qualifiers['protein_id'][0] = feature.qualifiers['protein_id'][0].replace(char, "_")

def _shorten_ids(idstring, options):
    contigstrmatch = re.search(r"onti?g?(\d+)\b", idstring)
    if not contigstrmatch:
        # if there is a substring "[Ss]caf(fold)XXX" use this number
        contigstrmatch = re.search(r"caff?o?l?d?(\d+)\b", idstring)
    if not contigstrmatch:
        # if there is a substring " cXXX" use this number
        contigstrmatch = re.search(r"\bc(\d+)\b", idstring)
    if contigstrmatch:
        contig_no = int(contigstrmatch.group(1))
    else:
        # if the contig number cannot be parsed out, just count the contigs from 1 to n
        contig_no = options.orig_record_idx

    return "c{ctg:05d}_{origid}..".format(ctg=contig_no, origid=idstring[:7])


def fix_record_name_id(seq_record, options):
    "Fix a seq record's name and id to be <= 16 characters, the GenBank limit; if record name is too long, add c000X prefix"

    if seq_record.id == "unknown.1":
        seq_record.id = "unk_seq_{ctg:05d}".format(ctg=options.orig_record_idx)
        logging.warn('Invalid sequence id "unknown.1", replaced by %s', seq_record.id)

    if seq_record.name == "unknown":
        seq_record.name = "unk_seq_{ctg:05d}".format(ctg=options.orig_record_idx)
        logging.warn('Invalid sequence name "unknown", replaced by %s', seq_record.name)

    if len(seq_record.id) > 16:
        oldid = seq_record.id

        #Check if it is a RefSeq accession number like NZ_AMZN01000079.1 that is just too long because of the version number behind the dot
        if (seq_record.id[-2] == "." and
                seq_record.id.count(".") == 1 and
                len(seq_record.id.partition(".")[0]) <= 16 and
                seq_record.id.partition(".")[0] not in options.all_record_ids):
            seq_record.id = seq_record.id.partition(".")[0]
            options.all_record_ids.add(seq_record.id)
        else: #Check if the ID suggested by _shorten_ids is unique
            if _shorten_ids(oldid, options) not in options.all_record_ids:
                seq_record.id = _shorten_ids(oldid, options)
                options.all_record_ids.add(seq_record.id)
            else:
                x = 0
                while "%s_%i" % (seq_record.id[:12], x) in options.all_record_ids:
                    x += 1
                seq_record.id = "%s_%i" % (seq_record.id[:12], x)
                options.all_record_ids.add(seq_record.id)

        logging.warn('Fasta header too long: renamed "%s" to "%s"', oldid, seq_record.id)
        if seq_record.id not in options.extrarecord:
            options.extrarecord[seq_record.id] = Namespace()
        if "extradata" not in options.extrarecord[seq_record.id]:
            options.extrarecord[seq_record.id].extradata = {}
        if "orig_id" not in options.extrarecord[seq_record.id].extradata:
            options.extrarecord[seq_record.id].extradata["orig_id"] = oldid
        #seq_record.id = "%s...%s" % (seq_record.id[:12], seq_record.id[-1])


    if len(seq_record.name) > 16:

        seq_record.name = _shorten_ids(seq_record.name, options)

        #seq_record.name = "%s...%s" % (seq_record.name[:12], seq_record.name[-1])

    if 'accession' in seq_record.annotations and \
       len(seq_record.annotations['accession']) > 16:
        acc = seq_record.annotations['accession']

        seq_record.annotations['accession'] = _shorten_ids(acc, options)
        # seq_record.annotations['accession'] = "%s...%s" % (acc[:12], acc[-1])

    # Remove illegal characters from name: otherwise, file cannot be written
    illegal_chars  = '''!"#$%&()*+,:;=>?@[]^`'{|}/ '''
    for char in seq_record.id:
        if char in illegal_chars:
            seq_record.id = seq_record.id.replace(char,"")
    for char in seq_record.name:
        if char in illegal_chars:
            seq_record.name = seq_record.name.replace(char,"")

def ascii_string(inputstring):
    return "".join([char for char in inputstring if char in (string.ascii_letters + string.digits + string.punctuation + string.whitespace)])

def get_gene_acc(feature):
    "Get the gene accesion number; if not defined, use content of locus_tag or gene qualifier instead"
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id'][0]

    # I needed to include this for clusterblast to work as well with the data from db (which all has protein_id) as with user data (which hasn't)
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    return "no_accession"

def get_gene_accs(feature):
    "Get array of the protein_id tags"
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id']
    elif 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag']
    elif 'gene' in feature.qualifiers:
        return feature.qualifiers['gene']
    return ["no_accession"]

def get_ncbi_gi(feature):
    """Get the NCBI gi from feature
    returns gi if found, None if not present"""

    gi = None
    if 'db_xref' in feature.qualifiers:
        for db_xref in feature.qualifiers['db_xref']:
            if 'GI:' in db_xref:
                gi = db_xref[3:]
    return gi

def get_gene_id(feature):
    "Get the gene ID from locus_tag, gene name or protein id, in that order"
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id'][0]
    return "no_tag_found"

def get_gene_annotation(feature):
    "Get the gene annotation from the product qualifier"
    if 'product' in feature.qualifiers:
        return feature.qualifiers['product'][0]
    return "unannotated orf"

def get_gene_accession(feature):
    "Get the gene accession number from protein_id, locus_tag or gene name, in that order"
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id'][0]
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    return "no_tag_found"

def get_full_path(current_file, file_to_add):
    "Get the full path of file_to_add in the same directory as current_file"
    return path.join(path.dirname(path.abspath(current_file)), file_to_add)

def get_feature_dict(seq_record):
    """Get a dictionary mapping features to their IDs"""
    features = get_cds_features(seq_record)
    feature_by_id = {}
    for feature in features:
        gene_id = get_gene_id(feature)
        feature_by_id[gene_id] = feature
    return feature_by_id

def get_feature_dict_protein_id(seq_record):
    """Get a dictionary mapping features to their protein_ids"""
    features = get_cds_features(seq_record)
    feature_by_id = {}
    for feature in features:
        gene_ids = get_gene_accs(feature)

        # As the HMMer result file does only contain result lines for 1 protein_id;
        # if an entry has 2 (or more) protein_id entries, lookups failed in the original version
        # Therefore I changed this method to use the protein_id list instead of just
        # extracting the [0] element - althouth we loose some memory here due to redundancy.
        for gene_id in gene_ids:
            feature_by_id[gene_id] = feature
    return feature_by_id

def get_multifasta(seq_record):
    """Extract multi-protein FASTA from all CDS features in sequence record"""
    features = get_cds_features(seq_record)
    all_fastas = []
    for feature in features:
        gene_id = get_gene_id(feature)
        fasta_seq = feature.qualifiers['translation'][0]
        if "-" in str(fasta_seq):
            fasta_seq = Seq(str(fasta_seq).replace("-",""), generic_protein)

        # Never write empty fasta entries
        if len(fasta_seq) == 0:
            logging.debug("No translation for %s, skipping", gene_id)
            continue

        all_fastas.append(">%s\n%s" % (gene_id, fasta_seq))
    full_fasta = "\n".join(all_fastas)
    return full_fasta

def get_specific_multifasta(features):
    """Extract multi-protein FASTA from provided features"""
    all_fastas = []
    for feature in features:
        gene_id = get_gene_id(feature)
        fasta_seq = feature.qualifiers['translation'][0]

        # Never write empty fasta entries
        if len(fasta_seq) == 0:
            logging.debug("No translation for %s, skipping", gene_id)
            continue

        all_fastas.append(">%s\n%s" % (gene_id, fasta_seq))
    full_fasta = "\n".join(all_fastas)
    return full_fasta

def get_aa_translation(seq_record, feature):
    """Obtain content for translation qualifier for specific CDS feature in sequence record"""
    fasta_seq = feature.extract(seq_record.seq).ungap('-').translate(to_stop=True)
    if len(fasta_seq) == 0:
        logging.debug("Retranslating %s with stop codons", feature.id)
        fasta_seq = feature.extract(seq_record.seq).ungap('-').translate()
    if "*" in str(fasta_seq):
        fasta_seq = Seq(str(fasta_seq).replace("*","X"), generic_protein)
    if "-" in str(fasta_seq):
        fasta_seq = Seq(str(fasta_seq).replace("-",""), generic_protein)

    return fasta_seq

def get_aa_sequence(feature, to_stop=False):
    """Extract sequence from specific CDS feature in sequence record"""
    fasta_seq = feature.qualifiers['translation'][0]
    if "*" in fasta_seq:
        if to_stop:
            fasta_seq = fasta_seq.split('*')[0]
        else:
            fasta_seq = fasta_seq.replace("*","X")
    if "-" in fasta_seq:
        fasta_seq = fasta_seq.replace("-","")
    return fasta_seq

def read_fasta(filename):
    # type: (str) -> Dict[str, str]
    """ reads a fasta file into a dict: id -> sequence, returns the dict """
    ids = []
    sequence_info = [] # type: List[str]
    with open(filename, "r") as fasta:
        current_seq = []
        for line in fasta:
            line = line.strip()
            if not line:
                continue
            if line[0] == '>':
                ids.append(line[1:].replace(" ", "_"))
                if current_seq:
                    sequence_info.append("".join(current_seq))
                    current_seq = []
            else:
                if not ids:
                    raise ValueError("Sequence before identifier in fasta file")
                if not line.replace("-","z").isalpha():
                    raise ValueError("Sequence contains non-alphabetic characters")
                current_seq.append(line)
    if current_seq:
        sequence_info.append("".join(current_seq))
    if len(ids) != len(sequence_info):
        raise ValueError("Fasta files contains different counts of sequences and ids")
    if not ids:
        logging.debug("Fasta file %s contains no sequences", filename)
        # TODO: refactor code using this func to deal with this ValueError instead
        # raise ValueError("Fasta file contains no sequences")
    return OrderedDict(zip(ids, sequence_info))

def writefasta(names, seqs, filename):
    "Write sequence to a file"
    e = 0
    f = len(names) - 1
    out_file = open(filename,"w")
    while e <= f:
        out_file.write(">")
        out_file.write(names[e])
        out_file.write("\n")
        out_file.write(seqs[e])
        out_file.write("\n")
        e += 1
    out_file.close()

def sortdictkeysbyvaluesrev(indict):
    items = [(value, key) for key, value in indict.items()]
    items.sort()
    items.reverse()
    return [key for value, key in items]

def get_genefinding_basedir(options):
    "Identify the basedir for glimmer files"
    basedir = ''
    if 'glimmer' in options:
        if 'basedir' in options.glimmer:
            basedir = options.glimmer.basedir
    return basedir


def get_overlaps_table(seq_record):
    """Given seq_record, returns array of overlapping genes
    and the corresponding gene_id-to-indexes"""
    overlaps = []
    overlap_by_id = {}
    features = get_cds_features(seq_record)
    if len(features) < 1:
        return overlaps, overlap_by_id
    features.sort(key=lambda feature: feature.location.start)
    i = 0
    j = i + 1
    cds_queue = []
    while i < len(features):
        if j >= len(features):
            break
        cds = features[i]
        ncds = features[j]
        if (cds.location.end <= ncds.location.start + 1):
            overlaps.append([])
            cds_queue.append(cds)
            for cds in cds_queue:
                overlap_by_id[get_gene_id(cds)] = len(overlaps) - 1
                overlaps[-1].append(cds)
            cds_queue = []
            i = j
        else:
            if (cds.location.end < ncds.location.end):
                cds_queue.append(cds)
                i = j
            else:
                cds_queue.append(ncds)
        j += 1
    overlaps.append([])
    cds_queue.append(features[i])
    for cds in cds_queue:
        overlap_by_id[get_gene_id(cds)] = len(overlaps) - 1
        overlaps[-1].append(cds)
    return overlaps, overlap_by_id


def zip_dir(dir_path, archive, prefix_to_remove=""):
    """Recursively add a directory's contents to a zip archive"""
    entries = os.listdir(dir_path)
    for entry in entries:
        entry = path.join(dir_path, entry)
        if path.isdir(entry):
            zip_dir(entry, archive, prefix_to_remove)
        else:
            arcname = entry.replace(prefix_to_remove + os.sep, "", 1)
            archive.write(entry, arcname)

def zip_path(dir_path, name):
    """Create a zip file for a given path"""

    # Don't zip in pre-existing zip files
    if path.exists(path.join(dir_path, name)):
        os.remove(path.join(dir_path, name))

    with TemporaryDirectory(change=True):
        try:
            archive = ZipFile(name, 'w', ZIP_DEFLATED)
            if path.isdir(dir_path):
                zip_dir(dir_path, archive, path.dirname(dir_path))
            else:
                archive.write(dir_path)
            archive.close()
        except LargeZipFile:
            archive = ZipFile(name, 'w', ZIP_DEFLATED, True)
            if path.isdir(dir_path):
                zip_dir(dir_path, archive, path.dirname(dir_path))
            else:
                archive.write(dir_path)
            archive.close()

        shutil.move(name, dir_path)

def log_status(status, level='running'):
    """Write status to the status file"""
    logging.info(status)

    options = get_config()
    if 'statusfile' not in options:
        return

    with open(options.statusfile, 'w') as statusfile:
        statusfile.write("%s: %s\n" % (level, status))


def get_git_version():
    """Get the sha1 of the current git version"""
    args = ['git', 'rev-parse', '--short', 'HEAD']
    try:
        out, _, _ = execute(args)
        return out.strip()
    except OSError:
        pass

    return ""


def get_version():
    """Get the current version string"""
    import antismash
    version = antismash.__version__
    git_version = get_git_version()
    if git_version != '':
        version += "-%s" % git_version

    return version

# this py2/py3 safe function will always cause a pylint error, so silence it
# pylint: disable=no-member
def create_protein_cleanup_table():
    from Bio.Data.IUPACData import protein_letters
    if hasattr(str, "maketrans"):
        maketrans = str.maketrans
    else:
        import string
        maketrans = string.maketrans
    return maketrans(protein_letters, ' ' * len(protein_letters))
# pylint: enable=no-member

def sort_XonY(x, y, reverse):
    """ Return a copy of x sorted according to the values in y
        :param x: the list to be sorted
        :param y: the list to use as sort values
        :param reverse: if True, sorts in decreasing order instead of increasing
        :return: the sorted list
    """
    return [i[0] for i in sorted(zip(x,y), key=lambda x: x[1], reverse=reverse)]


class RobustProteinAnalysis(ProteinAnalysis):
    TRANS_TABLE = create_protein_cleanup_table()

    def __init__(self, prot_sequence, monoisotopic=False, invalid="ignore"):
        if invalid not in ("ignore", "average"):
            raise ValueError("Invalid needs to be 'ignore' or 'average', not {!r}".format(invalid))
        self._invalid = invalid

        prot_sequence = prot_sequence.upper()

        self.original_sequence = prot_sequence
        prot_sequence = prot_sequence.translate(None, RobustProteinAnalysis.TRANS_TABLE)
        super(RobustProteinAnalysis, self).__init__(prot_sequence, monoisotopic)

    def molecular_weight(self):
        mw = super(RobustProteinAnalysis, self).molecular_weight()
        if self._invalid == "average":
            aa_difference = len(self.original_sequence) - len(self.sequence)
            mw += 110 * aa_difference

        return mw


def detect_unicode_in_file(filename):
    """"Detect unicode characters in an input file, raise """
    try:
        fh = codecs.open(filename, 'r', 'ascii')
        fh.read()
        fh.close()
    except ValueError:
        with codecs.open(filename, 'r', 'utf-8') as fh:
            lines = fh.readlines()
        for i, line in enumerate(lines):
            try:
                line.encode('ascii')
            except UnicodeEncodeError:
                return i, line
    return -1, None


def fetch_entries_from_ncbi(efetch_url):
    urltry = "n"
    nrtries = 0
    output = ""
    while urltry == "n" and nrtries < 4:
        try:
            nrtries += 1
            time.sleep(3)
            req = urllib2.Request(efetch_url)
            response = urllib2.urlopen(req)
            output = response.read()
            if len(output) > 5:
                urltry = "y"
        except (IOError, httplib.BadStatusLine, URLError, httplib.HTTPException):
            logging.error("Entry fetching from NCBI failed. Waiting for connection...")
            time.sleep(5)
    return output


def fix_wgs_master_record(seq_record):
    updated_seq_records = []
    # If seq_record is a WGS master record, parse out contig accession numbers and
    # download these as separate seq_records
    if 'wgs_scafld' in seq_record.annotations:
        contigranges = seq_record.annotations['wgs_scafld']
    else:
        contigranges = [seq_record.annotations['wgs']]
    allcontigs = []
    for contigrange in contigranges:
        if len(contigrange) == 1:
            allcontigs.extend(contigrange)
            continue
        startnumber, endnumber = '', ''
        alpha_tag = ''
        for char in contigrange[0].partition(".")[0]:
            if char.isdigit():
                startnumber += char
            else:
                alpha_tag += char
        for char in contigrange[1].partition(".")[0]:
            if char.isdigit():
                endnumber += char
        nrzeros = 0
        for char in startnumber:
            if char == "0":
                nrzeros += 1
            else:
                break
        contigrange = [alpha_tag + nrzeros * "0" + str(number) for number in xrange(int(startnumber), int(endnumber))]
        allcontigs.extend(contigrange)
    # Create contig groups of 50 (reasonable download size per download)
    nr_groups = len(allcontigs) / 50 + 1
    contig_groups = [list(np_array) for np_array in array_split(array(allcontigs), nr_groups)]
    # Return unchanged if no contigs supplied
    if len(contig_groups[0]) == 0 or (len(contig_groups[0]) == 1 and contig_groups[0][0] == ""):
        return [seq_record]
    # Download contigs and parse into seq_record objects
    for i, contig_group in enumerate(contig_groups):
        logging.info("Downloading contigs.. (%d-%d of %d)", 1+(i*50), min(50+(i*50), len(allcontigs)), len(allcontigs))
        efetch_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='
        efetch_url = efetch_url + ",".join(contig_group) + '&rettype=gbwithparts&retmode=text'
        output = fetch_entries_from_ncbi(efetch_url)
        if not len(output) > 5:
            break
        if "Resource temporarily unavailable" in output[:200] or "<h1>Server Error</h1>" in output[:500] \
                or "NCBI - WWW Error" in output[:500]:
            logging.error('ERROR: NCBI server temporarily unavailable: downloading contigs failed.')
            sys.exit(1)
        try:
            handle = StringIO(output)
            updated_seq_records.extend(list(SeqIO.parse(handle, 'genbank')))
        except ValueError as e:
            logging.error('Parsing %r failed: %s', "temporary contig file", e)
    return updated_seq_records


def fix_supercontig_record(seq_record):
    updated_seq_records = []
    # If seq_record is a supercontig record, reconstruct sequence and replace CONTIG feature by ORIGIN feature
    contig_info = seq_record.annotations['contig'].replace("join(", "")
    if contig_info[-1] == ")":
        contig_info = contig_info[:-1]
    allcontigparts = contig_info.split(",")
    accessions = []
    for part in contig_info.split(","):
        if "gap(" not in part:
            accessions.append(part.partition(":")[0].partition(".")[0].rpartition("complement(")[2].replace(")", ""))
    # Create contig groups of 50 (reasonable download size per download)
    nr_groups = len(accessions) / 50 + 1
    contig_groups = [list(np_array) for np_array in array_split(array(accessions), nr_groups)]
    # Return unchanged if no contigs supplied
    if len(contig_groups[0]) == 0 or (len(contig_groups[0]) == 1 and contig_groups[0][0] == ""):
        return [seq_record]
    # Download contig sequences based on their accessions
    contigseqdict = {}
    for i, contig_group in enumerate(contig_groups):
        logging.info("Downloading contigs.. (%d-%d of %d)", 1+(i*50), min(50+(i*50), len(allcontigparts)),
                     len(allcontigparts))
        efetch_url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id='
        efetch_url = efetch_url + ",".join(contig_group) + '&rettype=fasta&retmode=text'
        output = fetch_entries_from_ncbi(efetch_url)
        if not len(output) > 5:
            break
        if "Resource temporarily unavailable" in output[:200] or "<h1>Server Error</h1>" in output[:500]:
            logging.error('ERROR: NCBI server temporarily unavailable: downloading contigs failed.')
            sys.exit(1)
        sequences = [seq for seq in output.split(">") if len(seq) > 5]
        for sequence in sequences:
            for contig_acc in contig_group:
                if contig_acc in sequence.partition("\n")[0]:
                    contigseqdict[contig_acc] = sequence.partition("\n")[2].replace("\n", "")
    # Reconstruct supercontig sequence based on contig sequences and gaps
    fullsequence = ''
    for part in allcontigparts:
        if "gap(" in part:
            candidate_gap_int = part.partition('gap(')[2][:-1]
            if "unk" in candidate_gap_int:
                candidate_gap_int = candidate_gap_int.partition("unk")[2]
            if candidate_gap_int.isdigit():
                fullsequence += int(candidate_gap_int) * 'N'
            else:
                logging.error('Parsing supercontig file failed: faulty gap identifier' + part)
                sys.exit(1)
        else:
            accession = part.partition(":")[0].partition(".")[0].rpartition("complement(")[2]
            sequence = contigseqdict[accession]
            if "complement(" in part:
                sequence = str(Seq(sequence, generic_dna).reverse_complement())
            if ":" in part and ".." in part:
                seqstart, seqend = part.partition(":")[2].replace(")", "").replace("(", "").split("..")
                if int(seqstart) > 0:
                    seqstart = int(seqstart) - 1
                sequence = sequence[int(seqstart):int(seqend)]
            fullsequence += sequence
    # Add compiled sequence to seq_record
    seq = Seq(fullsequence, generic_dna)
    seq_record.seq = seq
    updated_seq_records.append(seq_record)
    return updated_seq_records


def process_wgs_master_scaffolds(seq_records):
    updated_seq_records = []
    n_wgs_scafld = 0
    n_supercontig = 0
    total_wgs_scafld = 0
    total_supercontig = 0
    for seq_record in seq_records:
        if 'wgs_scafld' in seq_record.annotations or 'wgs' in seq_record.annotations:
            total_wgs_scafld += 1
        elif 'contig' in seq_record.annotations:
            total_supercontig += 1
    for seq_record in seq_records:
        # We only need to download a sequence if our current seq is of type UnknownSeq
        run_download = isinstance(seq_record.seq, UnknownSeq)
        # Check if seq_record is a WGS master record or a supercontig record
        if run_download and ('wgs_scafld' in seq_record.annotations or 'wgs' in seq_record.annotations):
            n_wgs_scafld += 1
            logging.info("Fixing WGS master record, this might take a while.. "
                         "Internet connection is required. (%i/%i)",
                         n_wgs_scafld, total_wgs_scafld)
            updated_seq_records.extend(fix_wgs_master_record(seq_record))
        elif run_download and 'contig' in seq_record.annotations:
            n_supercontig += 1
            logging.info("Fixing Supercontig record, this might take a while.. "
                         "Internet connection is required. (%i/%i)",
                         n_supercontig, total_supercontig)
            updated_seq_records.extend(fix_supercontig_record(seq_record))
        else:
            if not run_download and ('wgs_scafld' in seq_record.annotations or
                                     'wgs' in seq_record.annotations or
                                     'contig' in seq_record.annotations):
                logging.info("Not downloading scaffolds for %s because sequence is present", seq_record.id)
            updated_seq_records.append(seq_record)

    return updated_seq_records
