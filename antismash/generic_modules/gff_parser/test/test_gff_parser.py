"Test suite for gff_parser module"

try:
    import unittest2
except ImportError:
    import unittest as unittest2
from os import path
from minimock import mock, restore, TraceTracker, assert_same_trace
from argparse import Namespace
from antismash import utils
from antismash import config
from antismash.generic_modules import gff_parser
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq


class FakeRecord(object):
    "class for generating a seq_record like data structure"
    def __init__(self, features=None, seq='AAAAAA'):
        if features is None:
            features = []
        self.features = features
        self.seq = Seq(seq)
    def __len__(self):
        return len(self.seq)


class FakeFeature(object):
    "class for generating a SeqFeature like datastructure"
    def __init__(self, type_, location=None, qualifiers=None):
        self.type = type_
        self.qualifiers = {} if qualifiers is None else qualifiers
        self.location = location

    def extract(self, seq):
        return seq

    def __repr__(self):
        return "FakeFeature(%r, %r, %r)" % (self.location, self.type,
                                            self.qualifiers)


class GffParserTest(unittest2.TestCase):
    def setUp(self):
        self.config = Namespace()
        config.set_config(self.config)
        self.config.gff3 = utils.get_full_path(__file__, "test_gff.gff")
        self.config.single_entries = False
        contig1 = FakeRecord(seq="".join(["A" for c in xrange(0, 2000)]))
        contig1.id = "CONTIG_1"
        contig2 = FakeRecord(seq="".join(["A" for c in xrange(0, 2000)]))
        contig2.id = "CONTIG_2"
        self.sequences = [contig1, contig2]

    def test_run(self):
        for sequence in self.sequences:
            gff_parser.run(sequence, self.config)
        len_cds_1 = len(utils.get_cds_features(self.sequences[0]))
        len_cds_2 = len(utils.get_cds_features(self.sequences[1]))
        detected_result = (len_cds_1, len_cds_2)
        expected_result = (1, 0)
        self.assertEqual(detected_result, expected_result, msg = "\nResult : %s\nExpected : %s" % (detected_result, expected_result))
