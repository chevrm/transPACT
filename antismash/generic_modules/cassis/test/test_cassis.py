# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2016 Thomas Wolf
# Leibniz Institute for Natural Product Research and
# Infection Biology -- Hans-Knoell-Institute (HKI)
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Test suite for the cassis cluster detection plugin"""

try:
    import unittest2
except ImportError:
    import unittest as unittest2

import os
import shutil
from argparse import Namespace
from shutil import copy
import re

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

from antismash import utils
from antismash.generic_modules import cassis


_this_dir = os.path.dirname(os.path.realpath(__file__))  # the folder of the file you are looking at


### helper methods ###
def convert_newline(string):
    """Convert all line endings to \n for OS independency"""
    if "\r\n" in string:
        return string.replace("\r\n", "\n") # win
    elif "\r" in string:
        return string.replace("\r", "\n") # mac
    else:
        return string


def read_generated_expected_file(generated_file, expected_file):
    """Read generated and expected files and save to string variables"""
    generated_string = ""
    with open(generated_file) as handle:
        generated_string = handle.read()
    generated_string = convert_newline(generated_string.rstrip())

    expected_string = ""
    with open(expected_file) as handle:
        expected_string = handle.read()
    expected_string = convert_newline(expected_string.rstrip())

    return [generated_string, expected_string]


def setup_fake_seq_record():
    """Set up a fake sequence record"""
    seq_record = SeqRecord(Seq(_fake_dna_seq)) # sequence --> see end of this file
    seq_record.name = "test"
    for i in range(1,10):
        feature = SeqFeature(type="gene")
        feature.qualifiers["locus_tag"] = ["gene" + str(i)]
        seq_record.features.append(feature)

    seq_record.features[0].location = FeatureLocation(100, 300, strand=1)
    seq_record.features[1].location = FeatureLocation(101, 299, strand=-1)
    seq_record.features[2].location = FeatureLocation(250, 350, strand=1)
    seq_record.features[3].location = FeatureLocation(500, 1000, strand=1)
    seq_record.features[4].location = FeatureLocation(1111, 1500, strand=-1)
    seq_record.features[5].location = FeatureLocation(2000, 2200, strand=-1)
    seq_record.features[6].location = FeatureLocation(2999, 4000, strand=1)
    seq_record.features[7].location = FeatureLocation(4321, 5678, strand=1)
    seq_record.features[8].location = FeatureLocation(6660, 9000, strand=-1)

    seq_record.features[3].qualifiers["sec_met"] = ["Type: faked"]
    seq_record.features[5].qualifiers["sec_met"] = ["Type: faked"]

    return seq_record


def setup_fake_options():
    """Set up a fake options namespace"""
    options = Namespace()
    options.outputfoldername = _this_dir
    options.cpus = 2
    options.skip_cleanup = False
    return options


### actual test classes and test methods ###
class TestCassisMethods(unittest2.TestCase):

    def tearDown(self):
        if os.path.exists(os.path.join(_this_dir, "meme")):
            shutil.rmtree(os.path.join(_this_dir, "meme"))
        if os.path.exists(os.path.join(_this_dir, "fimo")):
            shutil.rmtree(os.path.join(_this_dir, "fimo"))
        if os.path.isfile(os.path.join(_this_dir, "test_promoter_positions.csv")):
            os.remove(os.path.join(_this_dir, "test_promoter_positions.csv"))
        if os.path.isfile(os.path.join(_this_dir, "test_promoter_sequences.fasta")):
            os.remove(os.path.join(_this_dir, "test_promoter_sequences.fasta"))


    def test_plus_minus(self):
        self.assertEqual(len(cassis._plus_minus), 250)


    def test_get_anchor_genes(self):
        anchor_genes = ["gene4", "gene6"]
        seq_record = setup_fake_seq_record()
        self.assertEqual(cassis.get_anchor_genes(seq_record), anchor_genes)


    def test_ignore_overlapping(self):
        expected_not_ignored = ["gene1", "gene4", "gene5", "gene6", "gene7", "gene8", "gene9"]
        expected_ignored = ["gene2", "gene3"]
        seq_record = setup_fake_seq_record()
        not_ignored, ignored = cassis.ignore_overlapping(seq_record.features)

        self.assertEqual(map(lambda x: x.qualifiers["locus_tag"][0], ignored), expected_ignored)
        self.assertEqual(map(lambda x: x.qualifiers["locus_tag"][0], not_ignored), expected_not_ignored)


    def test_get_promoters(self):
        upstream_tss = 1000
        downstream_tss = 50
        options = setup_fake_options()
        seq_record = setup_fake_seq_record()
        genes, ignored_genes = cassis.ignore_overlapping(seq_record.features) # ignore ignored_genes

        # see cassis/promoterregions.png for details
        # [[start_prom1, end_prom1], [start_prom2, end_prom2], ...]
        expected_promoters = [
            [0, 150],
            [301, 550],
            [1450, 1999],
            [2150, 3049],
            [4001, 4371],
            [8950, 9603],
        ]

        promoters = cassis.get_promoters(seq_record, genes, upstream_tss, downstream_tss, options)
        self.assertEqual(map(lambda x: [x["start"], x["end"]], promoters), expected_promoters)

        # read expected files and save to string variable
        expected_sequences_file = ""
        with open(os.path.join(options.outputfoldername, "expected_promoter_sequences.fasta")) as handle:
            expected_sequences_file = handle.read()
        expected_sequences_file = convert_newline(expected_sequences_file.rstrip())

        expected_positions_file = ""
        with open(os.path.join(options.outputfoldername, "expected_promoter_positions.csv")) as handle:
            expected_positions_file = handle.read()
        expected_positions_file = convert_newline(expected_positions_file.rstrip())

        # read test files and save to string variable
        sequences_file = ""
        with open(os.path.join(options.outputfoldername, seq_record.name + "_promoter_sequences.fasta")) as handle:
            sequences_file = handle.read()
        sequences_file = convert_newline(sequences_file.rstrip())

        positions_file = ""
        with open(os.path.join(options.outputfoldername, seq_record.name + "_promoter_positions.csv")) as handle:
            positions_file = handle.read()
        positions_file = convert_newline(positions_file.rstrip())

        self.assertEqual(sequences_file, expected_sequences_file)
        self.assertEqual(positions_file, expected_positions_file)


    def test_get_anchor_promoter(self):
        anchor = "gene3"
        promoters = [
            {"id" : ["gene1"]},
            {"id" : ["gene2"]},
            {"id" : ["gene3", "gene4"]}, # sharing bidirectional promoter
            {"id" : ["gene5"]},
        ]
        self.assertEqual(cassis.get_anchor_promoter(anchor, promoters), 2)


    def test_get_promoter_sets(self):
        options = setup_fake_options()
        meme_dir = os.path.join(options.outputfoldername, "meme")
        anchor_promoter = 5
        promoters = [
            {"id" : ["gene1"],          "seq" : Seq("acgtacgtacgtacgt")},
            {"id" : ["gene2"],          "seq" : Seq("acgtacgtacgtacgt")},
            {"id" : ["gene3", "gene4"], "seq" : Seq("acgtacgtacgtacgt")},
            {"id" : ["gene5"],          "seq" : Seq("acgtacgtacgtacgt")},
            {"id" : ["gene6"],          "seq" : Seq("acgtacgtacgtacgt")},
            {"id" : ["gene7"],          "seq" : Seq("acgtacgtacgtacgt")}, # promoter with index=5 --> anchor promoter
            {"id" : ["gene8"],          "seq" : Seq("acgtacgtacgtacgt")},
            {"id" : ["gene9"],          "seq" : Seq("acgtacgtacgtacgt")},
        ]
        expected_promoter_sets = [
            {"plus" : 0, "minus" : 3, "score" : None},
            {"plus" : 0, "minus" : 4, "score" : None},
            {"plus" : 0, "minus" : 5, "score" : None},
            {"plus" : 1, "minus" : 2, "score" : None},
            {"plus" : 1, "minus" : 3, "score" : None},
            {"plus" : 1, "minus" : 4, "score" : None},
            {"plus" : 1, "minus" : 5, "score" : None},
            {"plus" : 2, "minus" : 1, "score" : None},
            {"plus" : 2, "minus" : 2, "score" : None},
            {"plus" : 2, "minus" : 3, "score" : None},
            {"plus" : 2, "minus" : 4, "score" : None},
            {"plus" : 2, "minus" : 5, "score" : None},
        ]
        self.assertEqual(cassis.get_promoter_sets(meme_dir, anchor_promoter, promoters), expected_promoter_sets)


    def test_filter_meme_results(self):
        options = setup_fake_options()
        meme_dir = os.path.join(options.outputfoldername, "meme")
        anchor = "AFUA_6G09660"
        promoter_sets = [
            {"plus" : 0, "minus" : 3},
        ]
        expected_motifs = [{
            "plus" : 0,
            "minus" : 3,
            "score" : "3.9e+003",
            "seqs" : [
                "TTTCGACCCGTC",
                "TTTCAAACCGTC",
                "TTTTGATTCGTC",
                "TTTTGACCGGTC",
                "TTTTAGACGGTC",
                "TTTTACCTCGTC",
                "TCTCGATCCGTC",
                "TTTCTATCCGTT",
                "TTTTGGACCGCC",
                "ATTTGGCCTGTC",
                "TGTTGTCTCGTC",
                "TTTGAGGCCGTC",
                "TTGTATTCTGTC",
                "TTTCTTCCTGTT",
            ]
        }]

        source = os.path.join(options.outputfoldername, "real_meme.xml") # this is a "real" MEME output file, I was too lazy to create my own fake XML file
        target = os.path.join(meme_dir, "+00_-03")
        if not os.path.exists(target):
            os.makedirs(target)
        copy(source, os.path.join(target, "meme.xml")) # overwrite meme.xml if exists

        self.assertEqual(cassis.filter_meme_results(meme_dir, promoter_sets, anchor), expected_motifs)
        binding_sites, expected_binding_sites = read_generated_expected_file(
            os.path.join(meme_dir, "+00_-03", "binding_sites.fasta"),
            os.path.join(options.outputfoldername, "expected_binding_sites.fasta")
        )
        self.assertEqual(binding_sites, expected_binding_sites)


    def test_filter_fimo_results(self):
        options = setup_fake_options() # necessary argument for the method but actually not used so far
        fimo_dir = os.path.join(options.outputfoldername, "fimo")
        motifs = [{"plus" : 0, "minus" : 3}]
        anchor_promoter = 1
        promoters = [
            {"start" : 10,  "end" : 14,  "id" : ["gene1"]},
            {"start" : 20,  "end" : 24,  "id" : ["gene2"]}, # anchor promoter
            {"start" : 30,  "end" : 34,  "id" : ["gene3"]},
            {"start" : 40,  "end" : 44,  "id" : ["gene4"]},
            {"start" : 50,  "end" : 54,  "id" : ["gene5"]},
            {"start" : 60,  "end" : 64,  "id" : ["gene6"]},
            {"start" : 70,  "end" : 74,  "id" : ["gene7"]},
            {"start" : 80,  "end" : 84,  "id" : ["gene8"]},
            {"start" : 90,  "end" : 94,  "id" : ["gene9"]},
            {"start" : 100, "end" : 104, "id" : ["gene10"]},
            {"start" : 110, "end" : 114, "id" : ["gene11"]},
            {"start" : 120, "end" : 124, "id" : ["gene12"]},
            {"start" : 130, "end" : 134, "id" : ["gene13"]},
            {"start" : 140, "end" : 144, "id" : ["gene14"]},
            {"start" : 150, "end" : 154, "id" : ["gene15"]},
            # need certain amount of promoters, otherwise the proportion of promoters with a motif (motif frequency) will be too high --> error
        ]
        expected_motifs = [{"plus" : 0, "minus" : 3, "hits" : {"gene1" : 1, "gene2" : 2}}]

        source = os.path.join(options.outputfoldername, "fake_short_fimo.txt") # fake FIMO output file, corresponding to expected_motifs
        target = os.path.join(fimo_dir, "+00_-03")
        if not os.path.exists(target):
            os.makedirs(target)
        copy(source, os.path.join(target, "fimo.txt")) # overwrite fimo.txt if exists

        self.assertEqual(cassis.filter_fimo_results(motifs, fimo_dir, promoters, anchor_promoter, options), expected_motifs)
        bs_per_promoter, expected_bs_per_promoter = read_generated_expected_file(
            os.path.join(target, "bs_per_promoter.csv"),
            os.path.join(options.outputfoldername, "expected_bs_per_promoter.csv")
        )
        self.assertEqual(bs_per_promoter, expected_bs_per_promoter)


    def test_get_islands(self):
        options = setup_fake_options
        motifs = [ # 2 different motifs
            {"plus" : 0, "minus" : 3, "hits" : {"gene1" : 1, "gene2" : 2}},
            {"plus" : 0, "minus" : 4, "hits" : {"gene2" : 3, "gene4" : 2, "gene5" : 1}},
        ]
        anchor_promoter = 1
        promoters = [
            {"start" : 10,  "end" : 14,  "id" : ["gene1"]},
            {"start" : 20,  "end" : 24,  "id" : ["gene2"]}, # anchor promoter
            {"start" : 30,  "end" : 34,  "id" : ["gene3"]},
            {"start" : 40,  "end" : 44,  "id" : ["gene4"]},
            {"start" : 50,  "end" : 54,  "id" : ["gene5"]},
            {"start" : 60,  "end" : 64,  "id" : ["gene6"]},
        ]
        expected_islands = [ # resulting in 2 different islands (this example)
            {
                # promoter (pos): 1 2 3 4 5 6
                # binding sites:  1 2 0 0 0 0
                # island:        |---|
                "start": {"start": 10, "end": 14, "id": ["gene1"]}, # first promoter of island
                "end":   {"start": 20, "end": 24, "id": ["gene2"]}, # last promoter of island
                "motif": {"hits": {"gene1": 1, "gene2": 2}, "plus": 0, "minus": 3}, # island is result of this motif
            },
            {
                # promoter (pos): 1 2 3 4 5 6
                # binding sites:  0 3 0 2 1 0
                # island:          |-------|
                "start": {"start": 20, "end": 24, "id": ["gene2"]}, # first promoter of island
                "end":   {"start": 50, "end": 54, "id": ["gene5"]}, # last promoter of island
                "motif": {"hits": {"gene2": 3, "gene4": 2, "gene5": 1}, "plus": 0, "minus": 4}, # island is result of this motif
            },
        ]
        self.assertEqual(cassis.get_islands(anchor_promoter, motifs, promoters, options), expected_islands)


    def test_sort_by_abundance(self):
        islands = [
            { # island 1: [gene1 -- gene2]
                "start": {"start": 1, "end": 1, "id": ["gene1"]},
                "end":   {"start": 2, "end": 2, "id": ["gene2"]},
                "motif": {"hits": {"gene1": 1, "gene2": 1},                                     "plus": 0, "minus": 3, "score" : 3},
            },
            { # island 2: [gene2 -- gene5]
                "start": {"start": 2, "end": 2, "id": ["gene2"]},
                "end":   {"start": 5, "end": 5, "id": ["gene5"]},
                "motif": {"hits": {"gene2": 1, "gene3": 1, "gene4": 1, "gene5": 1},             "plus": 3, "minus": 0, "score" : 2},
            },
            { # island 3: [gene1 -- gene5]
                "start": {"start": 1, "end": 1, "id": ["gene1"]},
                "end":   {"start": 5, "end": 5, "id": ["gene5"]},
                "motif": {"hits": {"gene1": 1, "gene2": 1, "gene3": 1, "gene4": 1, "gene5": 1}, "plus": 3, "minus": 3, "score" : 1},
            },
        ]
        # left border:  2x gene1, 1x gene2
        # right border: 2x gene5, 1x gene2

        expected_clusters = [
            { # cluster 1: [gene1 -- gene5] --> abundance 2+2 (most abundant)
                "start" : { "gene": "gene1", "abundance": 2, "score": 1, "plus": 3, "minus": 3 },
                "end":    { "gene": "gene5", "abundance": 2, "score": 1, "plus": 3, "minus": 3 },
            },
            { # cluster 3: [gene2 -- gene5] --> abundance 1+2, score 2+1 (better/lower)
                "start" : { "gene": "gene2", "abundance": 1, "score": 2, "plus": 3, "minus": 0 },
                "end":    { "gene": "gene5", "abundance": 2, "score": 1, "plus": 3, "minus": 3 },
            },
            { # cluster 2: [gene1 -- gene2] --> abundance 2+1, score 1+3 (worse, higher)
                "start" : { "gene": "gene1", "abundance": 2, "score": 1, "plus": 3, "minus": 3 },
                "end":    { "gene": "gene2", "abundance": 1, "score": 3, "plus": 0, "minus": 3 },
            },
            { # cluster 4: [gene2 -- gene2] --> abundance 1+1
                "start" : { "gene": "gene2", "abundance": 1, "score": 2, "plus": 3, "minus": 0 },
                "end":    { "gene": "gene2", "abundance": 1, "score": 3, "plus": 0, "minus": 3 },
            },
        ]
        # abundance: as high as possible
        # score: as low as possible

        self.assertEqual(cassis.sort_by_abundance(islands), expected_clusters)


    def test_check_cluster_predictions(self):
        seq_record = setup_fake_seq_record()
        promoters = [
            {"id" : ["gene1"]},
            {"id" : ["gene2"]},
            {"id" : ["gene3", "gene4"]},
        ]
        ignored_genes = [ # see captured logging
            SeqFeature(qualifiers = {"locus_tag" : ["gene2"]}),
        ]
        clusters = [
            { # only one cluster for testing
                "start" : { "gene": "gene1", "abundance": 1, "score": 1, "plus": 3, "minus": 3 },
                "end":    { "gene": "gene4", "abundance": 1, "score": 1, "plus": 3, "minus": 3 },
            },
        ]
        checked_clusters = [
            {
                "start" : { "gene": "gene1", "promoter" : "gene1",       "abundance": 1, "score": 1, "plus": 3, "minus": 3 },
                "end":    { "gene": "gene4", "promoter" : "gene3+gene4", "abundance": 1, "score": 1, "plus": 3, "minus": 3 },
                "genes":      4,
                "promoters" : 3,
            },
        ]

        self.assertEqual(cassis.check_cluster_predictions(clusters, seq_record, promoters, ignored_genes), checked_clusters)
        # -------------------- >> begin captured logging << --------------------
        # root: INFO: Best prediction (most abundant): 'gene1' -- 'gene4'
        # root: INFO: Upstream cluster border located at or next to sequence record border, prediction could have been truncated by record border
        # root: INFO: Downstream cluster border located at or next to sequence record border, prediction could have been truncated by record border
        # root: INFO: Gene 'gene2' is part of the predicted cluster, but it is overlapping with another gene and was ignored
        # root: INFO: Gene 'gene2' could have effected the cluster prediction
        # --------------------- >> end captured logging << ---------------------


    def test_cleanup_outdir(self):
        options = setup_fake_options()
        anchor_genes = ["gene1", "gene4"]
        cluster_predictions = {
            "gene1" : [
                { # only one anchor gene and cluster for testing
                    "start" : { "gene": "gene1", "promoter" : "gene1",       "abundance": 1, "score": 1, "plus": 3, "minus": 3 },
                    "end":    { "gene": "gene4", "promoter" : "gene3+gene4", "abundance": 1, "score": 1, "plus": 3, "minus": 3 },
                    "genes":      4,
                    "promoters" : 3,
                },
            ]
        }

        # create some empty test dirs, which should be deleted during the test
        os.makedirs(os.path.join(options.outputfoldername, "meme", "gene1", "+03_-03")) # prediction! --> keep!
        os.makedirs(os.path.join(options.outputfoldername, "fimo", "gene1", "+03_-03")) # prediction! --> keep!
        os.makedirs(os.path.join(options.outputfoldername, "meme", "gene1", "+04_-04")) # no prediction --> delete
        os.makedirs(os.path.join(options.outputfoldername, "fimo", "gene1", "+04_-04")) # no prediction --> delete
        os.makedirs(os.path.join(options.outputfoldername, "meme", "gene4", "+03_-03")) # no prediction --> delete
        os.makedirs(os.path.join(options.outputfoldername, "fimo", "gene4", "+03_-03")) # no prediction --> delete
        os.makedirs(os.path.join(options.outputfoldername, "meme", "gene4", "+04_-04")) # prediction for this gene, but not from this motif --> delete
        os.makedirs(os.path.join(options.outputfoldername, "fimo", "gene4", "+04_-04")) # prediction for this gene, but not from this motif --> delete

        cassis.cleanup_outdir(anchor_genes, cluster_predictions, options)

        # assert kept directories
        self.assertTrue("gene1" in os.listdir(os.path.join(options.outputfoldername, "meme")))
        self.assertTrue("gene1" in os.listdir(os.path.join(options.outputfoldername, "fimo")))
        self.assertTrue("+03_-03" in os.listdir(os.path.join(options.outputfoldername, "meme", "gene1")))
        self.assertTrue("+03_-03" in os.listdir(os.path.join(options.outputfoldername, "fimo", "gene1")))

        # assert deleted directories
        self.assertTrue("gene4" not in os.listdir(os.path.join(options.outputfoldername, "meme")))
        self.assertTrue("gene4" not in os.listdir(os.path.join(options.outputfoldername, "fimo")))
        self.assertTrue("+04_-04" not in os.listdir(os.path.join(options.outputfoldername, "meme", "gene1")))
        self.assertTrue("+04_-04" not in os.listdir(os.path.join(options.outputfoldername, "fimo", "gene1")))


class TestCassisHelperMethods(unittest2.TestCase):

    def test_mprint(self):
        plus = 3
        minus = 3
        self.assertEqual(cassis.mprint(3, 3), "+03_-03")


    def test_get_promoter_id(self):
        fake_promoter1 = {"id" : ["gene1"]}
        fake_promoter2 = {"id" : ["gene1", "gene2"]}
        self.assertEqual(cassis.get_promoter_id(fake_promoter1), "gene1")
        self.assertEqual(cassis.get_promoter_id(fake_promoter2), "gene1+gene2")


class TestCassisStorageMethods(unittest2.TestCase):

    def test_store_promoters(self):
        promoters = [
            {"id" : ["gene1"],          "seq" : Seq("cgtacgtacgt"), "start" : 10, "end" : 20},
            {"id" : ["gene2"],          "seq" : Seq("cgtacgtacgt"), "start" : 30, "end" : 40},
            {"id" : ["gene3", "gene4"], "seq" : Seq("cgtacgtacgt"), "start" : 50, "end" : 60},
        ]
        record_with_promoters = setup_fake_seq_record()
        cassis.store_promoters(promoters, record_with_promoters) # add ("store") promoters to seq_record

        record_without_promoters = setup_fake_seq_record() # just the same, without adding promoters

        # test if store_promoters changed any non-promoter feature (should not!)
        gene_features = len(record_without_promoters.features)
        for i in xrange(gene_features):
            self.assertEqual(repr(record_without_promoters.features[i]), repr(record_with_promoters.features[i]))

        # test promoter features
        self.assertEqual(len(record_without_promoters.features) + len(promoters), len(record_with_promoters.features))
        for i in xrange(len(promoters)):
            self.assertEqual(record_with_promoters.features[gene_features + i].type, "promoter")
            self.assertEqual(record_with_promoters.features[gene_features + i].qualifiers["seq"], ["cgtacgtacgt"])

        # expecially test bidirectional promoter feature (third promoter, last feature)
        self.assertEqual(record_with_promoters.features[gene_features + len(promoters) - 1].qualifiers["locus_tag"], ["gene3", "gene4"])
        self.assertEqual(record_with_promoters.features[gene_features + len(promoters) - 1].qualifiers["note"], ["bidirectional promoter"])


    def test_store_clusters(self):
        # this test is similar to test_store_promoters
        anchor = "gene3"
        clusters = [
            { # best prediction
                "start" : { "gene": "gene1", "promoter" : "gene1",       "abundance": 2, "score": 1, "plus": 3, "minus": 3 },
                "end":    { "gene": "gene4", "promoter" : "gene3+gene4", "abundance": 1, "score": 1, "plus": 3, "minus": 3 },
                "genes":      4,
                "promoters" : 3,
            },
            { # alternative prediction
                "start" : { "gene": "gene1", "promoter" : "gene1",       "abundance": 1, "score": 1, "plus": 4, "minus": 4 },
                "end":    { "gene": "gene5", "promoter" : "gene5",       "abundance": 1, "score": 1, "plus": 4, "minus": 4 },
                "genes":      4,
                "promoters" : 3,
            },
        ]
        record_with_clusters = setup_fake_seq_record()
        cassis.store_clusters(anchor, clusters, record_with_clusters)

        record_without_clusters = setup_fake_seq_record() # just the same, without adding clusters

        # test if store_clusters changed any non-cluster feature (should not!)
        gene_features = len(record_without_clusters.features)
        for i in xrange(gene_features):
            self.assertEqual(repr(record_without_clusters.features[i]), repr(record_with_clusters.features[i]))

        # test cluster features
        self.assertEqual(len(record_without_clusters.features) + len(clusters), len(record_with_clusters.features))
        for i in xrange(len(clusters)):
            self.assertEqual(record_with_clusters.features[gene_features + i].type, "cluster_border")
            self.assertEqual(record_with_clusters.features[gene_features + i].qualifiers["tool"], ["cassis"])
            self.assertEqual(record_with_clusters.features[gene_features + i].qualifiers["anchor"], [anchor])
            self.assertEqual(record_with_clusters.features[gene_features + i].qualifiers["genes"], [clusters[i]["genes"]])
            self.assertEqual(record_with_clusters.features[gene_features + i].qualifiers["promoters"], [clusters[i]["promoters"]])
            self.assertEqual(record_with_clusters.features[gene_features + i].qualifiers["gene_left"], [clusters[i]["start"]["gene"]])
            self.assertEqual(record_with_clusters.features[gene_features + i].qualifiers["gene_right"], [clusters[i]["end"]["gene"]])
            # don't test all feature qualifiers, only some


class TestCassisMainMethod(unittest2.TestCase):

    def tearDown(self):
        if os.path.exists(os.path.join(_this_dir, "meme")):
            shutil.rmtree(os.path.join(_this_dir, "meme"))
        if os.path.exists(os.path.join(_this_dir, "fimo")):
            shutil.rmtree(os.path.join(_this_dir, "fimo"))
        if os.path.isfile(os.path.join(_this_dir, "test_promoter_positions.csv")):
            os.remove(os.path.join(_this_dir, "test_promoter_positions.csv"))
        if os.path.isfile(os.path.join(_this_dir, "test_promoter_sequences.fasta")):
            os.remove(os.path.join(_this_dir, "test_promoter_sequences.fasta"))

    def test_detect(self):
        options = setup_fake_options()
        options.max_percentage = 100.0 # if lower, detect stops early
        seq_record = setup_fake_seq_record()
        cassis.detect(seq_record, options)
        # TODO it would be cool to capture the output (logging messages) of cassis.detect
        # one could than check if, e.g. "Best prediction (most abundant): 'gene1' -- 'gene9'" is in the output
        # but, if using nosetests, it cannot be captured, because nosetests already collects all output


class TestCassisUtils(unittest2.TestCase):

    def tearDown(self):
        if os.path.exists(os.path.join(_this_dir, "meme")):
            shutil.rmtree(os.path.join(_this_dir, "meme"))
        if os.path.exists(os.path.join(_this_dir, "fimo")):
            shutil.rmtree(os.path.join(_this_dir, "fimo"))
        if os.path.isfile(os.path.join(_this_dir, "test_promoter_sequences.fasta")):
            os.remove(os.path.join(_this_dir, "test_promoter_sequences.fasta"))

    def test_run_fimo(self):
        options = setup_fake_options()
        seq_record = setup_fake_seq_record()
        meme_dir = os.path.join(options.outputfoldername, "meme")
        fimo_dir = os.path.join(options.outputfoldername, "fimo")
        fimo_subdirs = ["+00_-03", "+03_-00"]

        # fimo input (i)
        #--> sequences to search in
        copy(
            os.path.join(options.outputfoldername, "expected_promoter_sequences.fasta"),
            os.path.join(options.outputfoldername, seq_record.name + "_promoter_sequences.fasta")
        )
        # fimo input (ii)
        # --> meme output (meme.html; motifs to search for)
        # --> file derived from meme output (binding_sites.fasta; only additional information for users; still checked by fimo)
        for path in fimo_subdirs:
            if not os.path.exists(os.path.join(meme_dir, path)):
                os.makedirs(os.path.join(meme_dir, path))
            copy(os.path.join(options.outputfoldername, "fake_meme.html"), os.path.join(meme_dir, path, "meme.html"))
            copy(os.path.join(options.outputfoldername, "fake_binding_sites.fasta"), os.path.join(meme_dir, path, "binding_sites.fasta"))

        utils.run_fimo(meme_dir, fimo_dir, seq_record, options)

        for path in fimo_subdirs:
            fimo_result, expected_fimo_result = read_generated_expected_file(
                                                    os.path.join(options.outputfoldername, "fimo", path, "fimo.txt"),
                                                    os.path.join(options.outputfoldername, "fake_long_fimo.txt"),
                                                )

            self.assertEqual(fimo_result, expected_fimo_result)


_fake_dna_seq = ("acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
    "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgtacgta"
)
