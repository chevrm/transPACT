from Bio.Seq import UnknownSeq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import unittest
# pylint: disable=unused-import
import helperlibs # used in mocking
# pylint: enable=unused-import
from minimock import mock, restore, TraceTracker, assert_same_trace, Mock
import antismash.output_modules.genbank as gb


def generate_cluster_feature(start, end, product, idx):
    """Generate cluster features for testing"""
    loc = FeatureLocation(start, end)
    notes = ["Cluster number: {}".format(idx)]
    cluster = SeqFeature(loc, type='cluster', qualifiers={'note': notes, 'product': [product]})
    return cluster





class TestGenbankWrite(unittest.TestCase):
    def setUp(self):
        self.tt = TraceTracker()
        mock('helperlibs.bio.seqio.write', tracker=self.tt)
        self.options = Mock('options', tracker=self.tt, outputfoldername='test', input_type='nucl')

        rec1_features = [
            generate_cluster_feature(23, 42, 'lantipeptide', 1),
            generate_cluster_feature(300, 500, 'nrps', 2)
        ]
        rec2_features = [
            generate_cluster_feature(50, 70, 'lassopeptide', 3),
            generate_cluster_feature(500, 700, 't1pks', 4)
        ]

        record1 = SeqRecord(UnknownSeq(1000), 'record1', name='record1', features=rec1_features)
        record2 = SeqRecord(UnknownSeq(1000), 'record2', name='record2', features=rec2_features)
        self.records = [record1, record2]

        self.expected_template = """\
Called helperlibs.bio.seqio.write(
    [SeqRecord(seq=UnknownSeq(19, alphabet = Alphabet(), character = '?'), id='record1', name='record1', description='<unknown description>', dbxrefs=[])],
    '{0}test/record1.cluster001.gbk',
    'genbank')
Called helperlibs.bio.seqio.write(
    [SeqRecord(seq=UnknownSeq(200, alphabet = Alphabet(), character = '?'), id='record1', name='record1', description='<unknown description>', dbxrefs=[])],
    '{0}test/record1.cluster002.gbk',
    'genbank')
Called helperlibs.bio.seqio.write(
    [SeqRecord(seq=UnknownSeq(20, alphabet = Alphabet(), character = '?'), id='record2', name='record2', description='<unknown description>', dbxrefs=[])],
    '{0}test/record1.cluster003.gbk',
    'genbank')
Called helperlibs.bio.seqio.write(
    [SeqRecord(seq=UnknownSeq(200, alphabet = Alphabet(), character = '?'), id='record2', name='record2', description='<unknown description>', dbxrefs=[])],
    '{0}test/record1.cluster004.gbk',
    'genbank')
Called helperlibs.bio.seqio.write(
    [SeqRecord(seq=UnknownSeq(1000, alphabet = Alphabet(), character = '?'), id='record1', name='record1', description='<unknown description>', dbxrefs=[]), SeqRecord(seq=UnknownSeq(1000, alphabet = Alphabet(), character = '?'), id='record2', name='record2', description='<unknown description>', dbxrefs=[])],
    '{0}test/record1.final.gbk',
    'genbank')"""

    def tearDown(self):
        restore()

    def test_write(self):
        "Test output_modules.genbank.write()"
        gb.write(self.records, self.options)
        expected = self.expected_template.format('')
        assert_same_trace(self.tt, expected)

    def test_write_abspath(self):
        "Test output_modules.genbank.write() with absolute path output folder name"
        self.options.outputfoldername = "/test"
        gb.write(self.records, self.options)
        expected = self.expected_template.format('/')
        assert_same_trace(self.tt, expected)
