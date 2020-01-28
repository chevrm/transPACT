"Test suite for utils.py library"

import unittest as unittest2
import shutil
import sys
import tempfile
import os
from os import path
from argparse import Namespace
from minimock import mock, restore, Mock, TraceTracker, assert_same_trace
from antismash import utils
from antismash import config
from Bio.SeqFeature import FeatureLocation
from Bio.Seq import UnknownSeq
# these are used in mocks
# pylint: disable=unused-import
import subprocess
import Bio.SearchIO
# pylint: enable=unused-import

# pylint doesn't like these, but I do
# pylint: disable=too-few-public-methods
class FakeRecord(object):
    "class for generating a seq_record like data structure"
    def __init__(self, features=None, seq='FAKESEQ', annotations=None):
        self.id = "FAKE"
        if features is None:
            features = []
        self.features = features
        self.seq = FakeSeq(seq)
        self.annotations = annotations if annotations else {}

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

class FakeSeq(object):
    "class for generating a Seq like datastructure"
    def __init__(self, seq):
        self.seq = seq

    def translate(self, dummy):
        return self.seq

    def __str__(self):
        return self.seq
# pylint: enable=too-few-public-methods


class TestUtils(unittest2.TestCase):

    def setUp(self):
        self.rec = FakeRecord()
        self.types = ["CDS", "gene", "CDS", "PFAM_domain", "cluster", "CDS", "aSDomain", ]
        idx = 1
        for t in self.types:
            f = FakeFeature(t)
            f.qualifiers['locus_tag'] = ["orf%04d" % idx]
            f.qualifiers['note'] = [
                'PFAM-Id: PF%05d' % idx,
                'smCOG: SMCOG{0:04d}:FAKE{0:04d} (Score: {0} ; E-value: 1.1e-1{0};'.format(idx)
            ]
            f.qualifiers['db_xref'] = ['PFAM: PF%05d' % idx]
            f.qualifiers['translation'] = ['FAKESEQ']
            self.rec.features.append(f)
            idx += 1
        self.trace_tracker = TraceTracker()

    def test_fasta_reader(self):
        def run_as_file(inputs):
            with tempfile.NamedTemporaryFile() as io_test:
                io_test.write("\n".join(inputs))
                io_test.flush()
                results = utils.read_fasta(io_test.name)
            return results
        # good lines
        res = run_as_file([">a","ABC",">b c","DEF"])
        self.assertEqual(len(res), 2)
        self.assertEqual(res["a"], "ABC")
        self.assertEqual(res["b_c"], "DEF")

        res = run_as_file([">a","ABC","D-F"])
        self.assertEqual(len(res), 1)
        self.assertEqual(res["a"], "ABCD-F")

        # bad lines
        for input_lines in [[">a"], [">a", "ABC", ">b"], ["DEF", ">a"], [">a", "7"]]:
            with self.assertRaises(ValueError):
                print(input_lines)
                run_as_file(input_lines)

    def test_get_all_features_of_type(self):
        "Test utils.get_all_features_of_type()"
        for t in set(self.types):
            f = utils.get_all_features_of_type(self.rec, t)
            self.assertEqual(len(f), self.types.count(t))

    def test_get_all_features_of_type_with_query(self):
        "Test utils.get_all_features_of_type_with_query()"
        f = utils.get_all_features_of_type_with_query(self.rec, 'aSDomain', 'note', 'PFAM-Id: PF00007')
        self.assertEqual(self.rec.features[6], f[0])

        f = utils.get_all_features_of_type_with_query(self.rec, 'aSDomain', 'note', 'PFAM-Id: PF00008')
        self.assertListEqual([], f)

    def test_get_cds_features(self):
        "Test utils.get_all_cds_features()"
        cds = utils.get_cds_features(self.rec)
        features = utils.get_all_features_of_type(self.rec, "CDS")
        self.assertListEqual(cds, features)


    def test_get_cluster_features(self):
        "Test utils.get_cluster_features()"
        clusters = utils.get_cluster_features(self.rec)
        features = utils.get_all_features_of_type(self.rec, "cluster")
        self.assertListEqual(clusters, features)


    def test_get_smcog_annotations(self):
        "Test utils.get_smcog_annotations()"
        expected_dict = {
            'orf0001': 'SMCOG0001',
            'orf0003': 'SMCOG0003',
            'orf0006': 'SMCOG0006',
        }
        expected_desc = {
            'SMCOG0001': 'FAKE0001 ',
            'SMCOG0003': 'FAKE0003 ',
            'SMCOG0006': 'FAKE0006 ',
        }
        smcogdict, smcog_desc = utils.get_smcog_annotations(self.rec)
        self.assertEqual(expected_dict, smcogdict)
        self.assertEqual(expected_desc, smcog_desc)


    def test_get_pfam_features(self):
        "Test utils.get_pfam_features()"
        motifs = utils.get_pfam_features(self.rec)
        features = utils.get_all_features_of_type(self.rec, "PFAM_domain")
        self.assertListEqual(motifs, features)

    def test_get_feature_dict(self):
        "Test utils.get_feature_dict()"
        fd = utils.get_feature_dict(self.rec)
        ids = [f.qualifiers['locus_tag'][0] for f in self.rec.features if f.type == "CDS"]
        keys = fd.keys()
        ids.sort()
        keys.sort()
        self.assertListEqual(ids, keys)

    def test_get_multifasta(self):
        "Test utils.get_multifasta"
        expected = """>orf0001
FAKESEQ
>orf0003
FAKESEQ
>orf0006
FAKESEQ"""

        ret = utils.get_multifasta(self.rec)
        self.assertMultiLineEqual(expected, ret)

    def test_get_overlaps_table(self):
        "Test utils.get_overlaps_table()"
        mock_features = [
            FakeFeature('CDS', FeatureLocation(10, 40), {"locus_tag" : ["G1"]}),
            FakeFeature('CDS', FeatureLocation(40, 50), {"locus_tag" : ["G2"]}),
            FakeFeature('CDS', FeatureLocation(45, 70), {"locus_tag" : ["G3"]}),
            FakeFeature('CDS', FeatureLocation(75, 100), {"locus_tag" : ["G4"]}),
            FakeFeature('CDS', FeatureLocation(101, 110), {"locus_tag" : ["G5"]}),
        ]
        mock_rec = FakeRecord(mock_features)
        result = utils.get_overlaps_table(mock_rec)
        expected = ([
                [mock_features[0]],
                [mock_features[1], mock_features[2]],
                [mock_features[3]],
                [mock_features[4]]
            ], {'G5': 3, 'G4': 2, 'G3': 1, 'G2': 1, 'G1': 0})
        self.assertEqual(result, expected, msg = result)


class TestContainmentHelpers(unittest2.TestCase):
    def setUp(self):
        self.features = [
            FakeFeature('cluster', FeatureLocation(25, 50), {'note': ['Cluster number: 1']}),
            FakeFeature('CDS', FeatureLocation(15, 20)),
            FakeFeature('PFAM_domain', FeatureLocation(15, 17)),
            FakeFeature('CDS', FeatureLocation(23, 42)),
            FakeFeature('CDS', FeatureLocation(45, 47)),
            FakeFeature('CDS', FeatureLocation(48, 55)),
            FakeFeature('aSDomain', FeatureLocation(4730, 4740)),
            FakeFeature('CDS', FeatureLocation(4700, 4710)),
            FakeFeature('CDS', FeatureLocation(4750, 4760)),
            FakeFeature('CDS', FeatureLocation(4790, 4812)),
            FakeFeature('cluster', FeatureLocation(4711, 4800), {'note': ['Cluster number: 2']}),
        ]

        self.record = FakeRecord(self.features)


    def test_get_cluster_cds_features(self):
        "Test utils.get_cluster_cds_features()"
        cluster1, cluster2 = utils.get_cluster_features(self.record)
        self.assertEqual(self.features[0], cluster1)
        self.assertEqual(self.features[-1], cluster2)

        clusterfeatures = utils.get_cluster_cds_features(cluster1, self.record)
        self.assertEqual(self.features[3:6], clusterfeatures)

        clusterfeatures = utils.get_cluster_cds_features(cluster2, self.record)
        self.assertEqual(self.features[-3:-1], clusterfeatures)


    def test_get_cluster_aSDomain_features(self):
        "Test utils.get_cluster_aSDomain_features()"
        cluster1, cluster2 = utils.get_cluster_features(self.record)
        self.assertEqual(self.features[0], cluster1)
        self.assertEqual(self.features[-1], cluster2)

        clusterfeatures = utils.get_cluster_aSDomain_features(cluster1, self.record)
        self.assertEqual([], clusterfeatures)

        clusterfeatures = utils.get_cluster_aSDomain_features(cluster2, self.record)
        self.assertEqual([self.features[-5]], clusterfeatures)


    def test_get_withincluster_cds_features(self):
        "Test utils.get_withincluster_cds_features()"
        features = utils.get_withincluster_cds_features(self.record)
        self.assertEqual(self.features[3:6] + self.features[-3:-1], features)


    def test_get_secmet_cds_featuers(self):
        "Test utils.get_secmet_cds_features()"
        self.features[3].qualifiers['sec_met'] = ["Type: Fake"]
        self.features[4].qualifiers['sec_met'] = ["Type: none"]
        features = utils.get_secmet_cds_features(self.record)
        self.assertEqual([self.features[3]], features)


    def test_get_pksnrps_cds_featuers(self):
        """Test utils.get_pksnrps_cds_features()"""
        self.features[3].qualifiers['sec_met'] = ["NRPS/PKS Domain: "]
        features = utils.get_pksnrps_cds_features(self.record)
        self.assertEqual([self.features[3]], features)


    def test_get_cluster_number(self):
        "Test utils.get_cluster_number()"
        # should return the actual number when it is present
        self.assertEqual(1, utils.get_cluster_number(self.features[0]))
        self.assertEqual(2, utils.get_cluster_number(self.features[-1]))

        # should return 0 otherwise
        no_number = FakeFeature('cluster', FeatureLocation(23, 42))
        self.assertEqual(0, utils.get_cluster_number(no_number))


    def test_get_sorted_cluster_features(self):
        "Test utils.get_sorted_cluster_features()"
        res = utils.get_sorted_cluster_features(self.record)
        self.assertEqual([self.features[0], self.features[-1]], res)

        # remove both cluster records
        self.features.pop(0)
        self.features.pop()
        self.assertEqual([], utils.get_sorted_cluster_features(self.record))


    def test_get_cluster_type(self):
        "Test utils.get_cluster_type()"
        cluster = FakeFeature('cluster', FeatureLocation(23, 42),
                              {'product': ['fake']})
        self.assertEqual('fake', utils.get_cluster_type(cluster))


    def test_get_structure_pred(self):
        "Test utils.get_structure_pred()"
        cluster = FakeFeature('cluster', FeatureLocation(23, 42),
                              {'product': ['fake']})

        self.assertEqual('N/A', utils.get_structure_pred(cluster))

        cluster.qualifiers['product'][0] = 'ectoine'
        self.assertEqual('ectoine', utils.get_structure_pred(cluster))

        cluster.qualifiers['note'] = ['Monomers prediction: fake']
        self.assertEqual('fake', utils.get_structure_pred(cluster))


    def test_get_cluster_features_of_type(self):
        "Test utils.get_cluster_features_of_type()"
        self.features[0].qualifiers['product'] = ['nrps']
        self.features[-1].qualifiers['product'] = ['lanthipeptide']

        self.assertEqual([], utils.get_cluster_features_of_type(self.record, 'pks'))
        self.assertEqual([self.features[0]],
                         utils.get_cluster_features_of_type(self.record, 'nrps'))
        self.assertEqual([self.features[-1]],
                         utils.get_cluster_features_of_type(self.record, 'lanthipeptide'))


    def test_features_overlap(self):
        "Test utils.features_overlap()"
        self.assertFalse(utils.features_overlap(self.features[0], self.features[1]))
        self.assertFalse(utils.features_overlap(self.features[1], self.features[0]))

        self.assertTrue(utils.features_overlap(self.features[0], self.features[3]))
        self.assertTrue(utils.features_overlap(self.features[3], self.features[0]))

        self.assertTrue(utils.features_overlap(self.features[0], self.features[4]))
        self.assertTrue(utils.features_overlap(self.features[4], self.features[0]))

        self.assertTrue(utils.features_overlap(self.features[0], self.features[5]))
        self.assertTrue(utils.features_overlap(self.features[5], self.features[0]))


    def check_sorting(self):
        for i in range(len(self.features) - 1):
            if self.features[i].location.start > \
               self.features[i+1].location.start:
                return False
            if self.features[i].location.start == \
               self.features[i+1].location.start and \
               self.features[i].location.end > \
               self.features[i+1].location.end:
                return False

        return True

    def test_sort_features(self):
        "Test utils.sort_features()"

        self.assertFalse(self.check_sorting())

        utils.sort_features(self.record)

        self.assertTrue(self.check_sorting())


    def test_cmp_feature_location(self):
        "Test utils.cmp_feature_location()"
        self.assertEqual(0, utils.cmp_feature_location(self.features[0], self.features[0]))
        self.assertEqual(-1, utils.cmp_feature_location(self.features[1], self.features[0]))
        self.assertEqual(1, utils.cmp_feature_location(self.features[0], self.features[1]))


class TestExecutableHelpers(unittest2.TestCase):
    def setUp(self):
        self.config = Namespace()
        self.config.cpus = 2
        config.set_config(self.config)
        self.tt = TraceTracker()
        proc = Mock('proc', tracker=self.tt, returncode=0)
        proc.communicate = Mock('proc.communicate', returns=('output', 'error'), tracker=self.tt)
        mock('subprocess.Popen', tracker=self.tt, returns=proc)
        self.tmpdir = tempfile.mkdtemp(prefix="as_tests_util")

    def tearDown(self):
        restore()
        shutil.rmtree(self.tmpdir)

    def test_locate_executable(self):
        "Test utils.locate_executable()"
        test_executable = sys.argv[0]
        if not path.isfile(test_executable):
            self.skipTest("%r is not a file, skipping test" % test_executable)

        if not os.access(test_executable, os.X_OK):
            self.skipTest("%r is not executable, skipping test" % test_executable)

        self.assertEqual(utils.locate_executable(test_executable),
                         test_executable)

        short_path = path.split(test_executable)[1]
        self.assertEqual(utils.locate_executable(short_path), test_executable)

        self.assertIsNone(utils.locate_executable("totallymaedup"))


    def test_locate_file(self):
        "Test utils.locate_file()"
        self.assertEqual(__file__, utils.locate_file(__file__))
        self.assertIsNone(utils.locate_file('totallymadeup'))


    def test_excute_no_input(self):
        "Test utils.execute() without stdin input"
        expected = """    Called subprocess.Popen(
        ['fake', '--with', 'parameters'],
        stderr=-1,
        stdin=None,
        stdout=-1)
    Called proc.communicate(input=None)"""
        cmd = ['fake', '--with', 'parameters']

        out, err, retcode = utils.execute(cmd)

        self.assertEqual('output', out)
        self.assertEqual('error', err)
        self.assertEqual(0, retcode)
        assert_same_trace(self.tt, expected)


    def test_execute_with_input(self):
        "Test utils.execute() with stdin input"
        expected = """    Called subprocess.Popen(
        ['fake', '--with', 'parameters'],
        stderr=-1,
        stdin=-1,
        stdout=-1)
    Called proc.communicate(input='fake input')"""
        cmd = ['fake', '--with', 'parameters']

        out, err, retcode = utils.execute(cmd, input='fake input')

        self.assertEqual('output', out)
        self.assertEqual('error', err)
        self.assertEqual(0, retcode)
        assert_same_trace(self.tt, expected)


    def test_run_hmmsearch(self):
        "Test utils.run_hmmsearch()"
        mock('Bio.SearchIO.parse', tracker=self.tt, returns=['mock result'])
        mock('utils.execute', tracker=self.tt, returns=('output', 'error', 0))

        expected = r"""    Called utils.execute(
        ['hmmsearch', '--cpu', '2', '-o', '/dev/null', '--domtblout', 'result.domtab', 'fake.hmm', '-'],
        input='>testinput\nMADEUP')
    Called Bio.SearchIO.parse(
        'result.domtab', 'hmmsearch3-domtab')"""

        hits = utils.run_hmmsearch('fake.hmm', ">testinput\nMADEUP")
        self.assertEqual(len(hits), 1)

        hit = hits.pop()
        self.assertEqual('mock result', hit)
        assert_same_trace(self.tt, expected)


    def test_run_hmmsearch_nothreads(self):
        "Test utils.run_hmmsearch() without multithreading"
        mock('Bio.SearchIO.parse', tracker=self.tt, returns=['mock result'])
        mock('utils.execute', tracker=self.tt, returns=('output', 'error', 0))
        self.config.hmmer3 = Namespace(multithreading=False)

        expected = r"""    Called utils.execute(
        ['hmmsearch', '-o', '/dev/null', '--domtblout', 'result.domtab', 'fake.hmm', '-'],
        input='>testinput\nMADEUP')
    Called Bio.SearchIO.parse(
        'result.domtab', 'hmmsearch3-domtab')"""

        hits = utils.run_hmmsearch('fake.hmm', ">testinput\nMADEUP")
        self.assertEqual(len(hits), 1)

        hit = hits.pop()
        self.assertEqual('mock result', hit)
        assert_same_trace(self.tt, expected)


    def test_run_hmmscan(self):
        "Test utils.run_hmmscan()"
        mock('Bio.SearchIO.parse', tracker=self.tt, returns=['mock result'])
        mock('utils.execute', tracker=self.tt, returns=('output', 'error', 0))

        expected = r"""    Called utils.execute(
        ['hmmscan', '--cpu', '2', '--nobias', 'fake.hmm', '-'],
        input='>testinput\nMADEUP')
    Called Bio.SearchIO.parse(
        <cStringIO.StringI object at ...>,
        'hmmer3-text')"""

        hits = utils.run_hmmscan('fake.hmm', ">testinput\nMADEUP")
        self.assertEqual(len(hits), 1)

        hit = hits.pop()
        self.assertEqual('mock result', hit)
        assert_same_trace(self.tt, expected)


    def test_run_hmmscan_nothreads(self):
        "Test utils.run_hmmscan() without multithreading"
        mock('Bio.SearchIO.parse', tracker=self.tt, returns=['mock result'])
        mock('utils.execute', tracker=self.tt, returns=('output', 'error', 0))
        self.config.hmmer3 = Namespace(multithreading=False)

        expected = r"""    Called utils.execute(
           ['hmmscan', '--nobias', 'fake.hmm', '-'],
           input='>testinput\nMADEUP')
    Called Bio.SearchIO.parse(
           <cStringIO.StringI object at ...>,
           'hmmer3-text')"""

        hits = utils.run_hmmscan('fake.hmm', ">testinput\nMADEUP")
        self.assertEqual(len(hits), 1)

        hit = hits.pop()
        self.assertEqual('mock result', hit)
        assert_same_trace(self.tt, expected)

    def test_run_hmmscan_write_resultfile(self):
        """Test utils.run_hmmscan() writing a results file"""
        mock('Bio.SearchIO.parse', tracker=self.tt, returns=['mock result'])
        mock('utils.execute', tracker=self.tt, returns=('output', 'error', 0))

        expected = r"""    Called utils.execute(
        ['hmmscan', '--cpu', '2', '--nobias', 'fake.hmm', '-'],
        input='>testinput\nMADEUP')
    Called Bio.SearchIO.parse(
        <cStringIO.StringI object at ...>,
        'hmmer3-text')"""

        results_file = path.join(self.tmpdir, 'fake_hmmscan_output.txt')
        hits = utils.run_hmmscan('fake.hmm', ">testinput\nMADEUP", results_file=results_file)
        self.assertEqual(len(hits), 1)

        hit = hits.pop()
        self.assertEqual('mock result', hit)
        assert_same_trace(self.tt, expected)
        self.assertTrue(path.exists(results_file))
        self.assertEqual(open(results_file).read(), 'output')

    def test_run_hmmpfam2(self):
        "Test utils.run_hmmpfam2()"
        mock('Bio.SearchIO.parse', tracker=self.tt, returns=['mock result'])
        mock('utils.execute', tracker=self.tt, returns=('output', 'error', 0))

        expected = r"""    Called utils.execute(
        ['hmmpfam2', '--cpu', '2', 'fake.hmm', '-'],
        input='>testinput\nMADEUP')
    Called Bio.SearchIO.parse(
        <cStringIO.StringI object at ...>,
        'hmmer2-text')"""

        hits = utils.run_hmmpfam2('fake.hmm', ">testinput\nMADEUP")
        self.assertEqual(len(hits), 1)

        hit = hits.pop()
        self.assertEqual('mock result', hit)
        assert_same_trace(self.tt, expected)


    def test_run_hmmpfam2_nothreads(self):
        "Test utils.run_hmmpfam2() without multithreading"
        mock('Bio.SearchIO.parse', tracker=self.tt, returns=['mock result'])
        mock('utils.execute', tracker=self.tt, returns=('output', 'error', 0))
        self.config.hmmer2 = Namespace(multithreading=False)

        expected = r"""    Called utils.execute(
        ['hmmpfam2', 'fake.hmm', '-'],
        input='>testinput\nMADEUP')
    Called Bio.SearchIO.parse(
        <cStringIO.StringI object at ...>,
        'hmmer2-text')"""

        hits = utils.run_hmmpfam2('fake.hmm', ">testinput\nMADEUP")
        self.assertEqual(len(hits), 1)

        hit = hits.pop()
        self.assertEqual('mock result', hit)
        assert_same_trace(self.tt, expected)


    def test_get_full_path(self):
        "Test utils.get_full_path()"
        expected = path.join(path.dirname(path.abspath(__file__)), 'fake')
        ret = utils.get_full_path(__file__, 'fake')

        self.assertEqual(ret, expected)


    def test_get_git_vesion(self):
        "Test utils.get_git_version()"
        mock('utils.execute', tracker=self.tt, returns=('deadbeef', '', 0))

        expected_trace = r"""  Called utils.execute(['git', 'rev-parse', '--short', 'HEAD'])
  Called utils.execute(['git', 'rev-parse', '--short', 'HEAD'])"""

        ret = utils.get_git_version()
        self.assertEqual(ret, 'deadbeef')

        mock('utils.execute', tracker=self.tt, raises=OSError)
        ret = utils.get_git_version()
        self.assertEqual(ret, '')

        assert_same_trace(self.tt, expected_trace)


class TestAnalysisHelpers(unittest2.TestCase):
    def test_get_gene_id_locus_tag(self):
        "Test utils.get_gene_id() with locus tag"
        expected = 'test_tag'
        f = FakeFeature("CDS")
        f.qualifiers['locus_tag'] = [expected]

        ret = utils.get_gene_id(f)
        self.assertEqual(ret, expected)

    def test_get_ncbi_gi(self):
        "Test utils.get_ncbi_gi() with gi tag"

        f = FakeFeature('CDS')
        f.qualifiers['db_xref'] = ['GI:test_gi']

        ret = utils.get_ncbi_gi(f)
        self.assertEqual(ret, 'test_gi')

        del f.qualifiers['db_xref']
        ret = utils.get_ncbi_gi(f)
        self.assertEqual(ret, None)

    def test_get_gene_id_gene(self):
        "Test utils.get_gene_id() with gene tag"
        expected = 'test_gene'
        f = FakeFeature("CDS")
        f.qualifiers['gene'] = [expected]

        ret = utils.get_gene_id(f)
        self.assertEqual(ret, expected)



    def test_get_gene_id_protein_id(self):
        "Test utils.get_gene_id() with protein_id tag"
        expected = 'test_id'
        f = FakeFeature("CDS")
        f.qualifiers['protein_id'] = [expected]

        ret = utils.get_gene_id(f)
        self.assertEqual(ret, expected)

    def test_get_gene_id_locus_tag_gene(self):
        "Test utils.get_gene_id() with locus tag and gene tag"
        expected = 'test_tag'
        f = FakeFeature("CDS")
        f.qualifiers['locus_tag'] = [expected]
        f.qualifiers['gene'] = ['test_gene']

        ret = utils.get_gene_id(f)
        self.assertEqual(ret, expected)

    def test_get_gene_id_locus_tag_protein(self):
        "Test utils.get_gene_id() with locus tag and protein_id tag"
        expected = 'test_tag'
        f = FakeFeature("CDS")
        f.qualifiers['locus_tag'] = [expected]
        f.qualifiers['protein_id'] = ['test_id']

        ret = utils.get_gene_id(f)
        self.assertEqual(ret, expected)

    def test_get_gene_id_gene_protein(self):
        "Test utils.get_gene_id() with gene tag and protein_id tag"
        expected = 'test_gene'
        f = FakeFeature("CDS")
        f.qualifiers['gene'] = [expected]
        f.qualifiers['protein_id'] = ['test_id']

        ret = utils.get_gene_id(f)
        self.assertEqual(ret, expected)

    def test_get_gene_id_no_id(self):
        "Test utils.get_gene_id() without any useable id"
        expected = 'no_tag_found'
        f = FakeFeature("CDS")

        ret = utils.get_gene_id(f)
        self.assertEqual(ret, expected)


    def test_get_gene_annotation(self):
        "Test utils.get_gene_annotation()"
        feature = FakeFeature("CDS")
        self.assertEqual('unannotated orf', utils.get_gene_annotation(feature))

        feature.qualifiers['product'] = ['fake']
        self.assertEqual('fake', utils.get_gene_annotation(feature))


    def test_get_gene_accession_protein(self):
        "Test utils.get_gene_accession() with protein_id"
        feature = FakeFeature("CDS")
        feature.qualifiers['protein_id'] = ['fake007']
        self.assertEqual('fake007', utils.get_gene_accession(feature))


    def test_get_gene_accession_locus(self):
        "Test utils.get_gene_accession() with locus_tag"
        feature = FakeFeature("CDS")
        feature.qualifiers['locus_tag'] = ['fake007']
        self.assertEqual('fake007', utils.get_gene_accession(feature))


    def test_get_gene_accession_gene(self):
        "Test utils.get_gene_accession() with gene"
        feature = FakeFeature("CDS")
        feature.qualifiers['gene'] = ['fake007']
        self.assertEqual('fake007', utils.get_gene_accession(feature))


    def test_get_gene_accession_no_id(self):
        "Test utils.get_gene_accession() without available id"
        feature = FakeFeature("CDS")
        self.assertEqual('no_tag_found', utils.get_gene_accession(feature))


    def test_get_aa_sequence(self):
        "Test utils.get_aa_sequence() for straightforward translation"
        expected = 'MAGIC'
        f = FakeFeature("CDS")
        f.qualifiers['translation'] = [expected]

        ret = utils.get_aa_sequence(f)
        self.assertEqual(expected, ret)


    def test_get_aa_sequence_stop(self):
        "Test utils.get_aa_sequence() for translation including a stop codon"
        inseq = 'MAGIC*'
        expected = 'MAGICX'
        f = FakeFeature("CDS")
        f.qualifiers['translation'] = [inseq]

        ret = utils.get_aa_sequence(f)
        self.assertEqual(expected, ret)


    def test_get_aa_sequence_gap(self):
        "Test utils.get_aa_sequence() for translation including a gap"
        inseq = 'MA-GIC'
        expected = 'MAGIC'
        f = FakeFeature("CDS")
        f.qualifiers['translation'] = [inseq]

        ret = utils.get_aa_sequence(f)
        self.assertEqual(expected, ret)


    def test_get_aa_sequence_to_stop(self):
        "Test utils.get_aa_sequence() for translation up to a stop codon"
        inseq = 'MAGIC*SEQ'
        expected = 'MAGIC'
        f = FakeFeature("CDS")
        f.qualifiers['translation'] = [inseq]

        ret = utils.get_aa_sequence(f, to_stop=True)
        self.assertEqual(expected, ret)

    def test_fix_record_name_id(self):
        """Test utils.fix_record_name_id"""
        record = FakeRecord()
        options = Namespace()
        record.id = "GoodID.1"
        record.name = "GoodID"
        record.annotations = {}
        options.all_record_ids = {"GoodID.1"}
        options.orig_record_idx = 1
        options.extrarecord = {}

        utils.fix_record_name_id(record, options)
        self.assertEqual(record.id, "GoodID.1")

        # No annotation at all results in the 'unknown' name/id
        record.name = "unknown"
        record.id = "unknown.1"
        options.all_record_ids = set()
        utils.fix_record_name_id(record, options)
        self.assertEqual(record.id, "unk_seq_00001")
        self.assertEqual(record.name, "unk_seq_00001")

        # Too long ID check for valid WGS IDs
        new_name = "NZ_AMZN01000079"
        record.id = "NZ_AMZN01000079.1"
        record.name = "NZ_AMZN01000079"
        options.all_record_ids = set()
        options.extrarecord = {}
        utils.fix_record_name_id(record, options)
        self.assertEqual(record.id, new_name)
        self.assertEqual(record.name, new_name)
        self.assertTrue(new_name in options.all_record_ids)
        self.assertTrue(new_name in options.extrarecord)

        # Too long ID for things that are not valid WGS IDs
        new_name = "c12345_MySuper.."
        record.id = "MySuperLongContig0000012345"
        record.name = "MySuperLongContig0000012345"
        options.all_record_ids = set()
        options.extrarecord = {}
        utils.fix_record_name_id(record, options)
        self.assertEqual(record.id, new_name)
        self.assertEqual(record.name, new_name)
        self.assertTrue(new_name in options.all_record_ids)
        self.assertTrue(new_name in options.extrarecord)

        # Too long ID with a conflict
        new_name = "MySuperLongC_1"
        record.id = "MySuperLongContig0000012345"
        record.name = "MySuperLongContig0000012345"
        options.all_record_ids = {"c12345_MySuper..", "MySuperLongC_0"}
        options.extrarecord = {}
        utils.fix_record_name_id(record, options)
        self.assertEqual(record.id, new_name)
        # record name not unique?
        self.assertEqual(record.name, "c12345_MySuper..")
        self.assertTrue(new_name in options.all_record_ids)
        self.assertTrue(new_name in options.extrarecord)


class TestRobustProteinAnalysis(unittest2.TestCase):

    def test_create_protein_cleanup_table(self):
        """Test utils.create_protein_cleanup_table() is sane"""
        trans_table = utils.create_protein_cleanup_table()
        orig_string = "MAGICXHAT"
        expected_deletions = "     X   "
        expected_result = "MAGICHAT"
        deletions = orig_string.translate(trans_table)
        result = orig_string.translate(None, trans_table)
        self.assertEqual(expected_deletions, deletions)
        self.assertEqual(expected_result, result)

    def test_init(self):
        """Test RobustProteinAnalysis initialisation"""
        rpa = utils.RobustProteinAnalysis("MAGICHAT", invalid='ignore')
        self.assertIsInstance(rpa, utils.RobustProteinAnalysis)

        rpa = utils.RobustProteinAnalysis("MAGICHAT", invalid='average')
        self.assertIsInstance(rpa, utils.RobustProteinAnalysis)

        self.assertRaises(ValueError, utils.RobustProteinAnalysis, "MAGICHAT", invalid='funky')

    def test_uppercase(self):
        """Test RobustProteinAnalysis converts passed sequence to upper case"""
        rpa = utils.RobustProteinAnalysis("Magichat")
        self.assertEqual("MAGICHAT", rpa.original_sequence)
        self.assertEqual("MAGICHAT", str(rpa.sequence))

    def test_molecular_weight_ignore(self):
        """Test RobustProteinAnalysis.molecular_weight() calculates correct weight with invalid='ignore'"""
        rpa = utils.RobustProteinAnalysis("MAGICXHAT")
        self.assertEqual(802.9621, rpa.molecular_weight())

    def test_molecular_weight_average(self):
        """Test RobustProteinAnalysis.molecular_weight() calculates correct weight with invalid='average'"""
        rpa = utils.RobustProteinAnalysis("MAGICXHAT", invalid='average')
        self.assertEqual(912.9621, rpa.molecular_weight())


class TestScaffoldDownload(unittest2.TestCase):
    """Test WGS/WGS_SCAFLD/CONTIG"""

    def setUp(self):
        self.tt = TraceTracker()
        # Fall over when we try to call the real download functions
        mock('utils.fix_supercontig_record', tracker=self.tt, raises=NotImplementedError)
        mock('utils.fix_wgs_master_record', tracker=self.tt, raises=NotImplementedError)

    def tearDown(self):
        restore()

    def test_process_wgs_master_scaffolds_no_scaffolds(self):
        """Test process_wgs_master_scaffolds() does not change the list of records if there are no scaffolds"""
        fake_records = [FakeRecord(), FakeRecord()]
        procesed = utils.process_wgs_master_scaffolds(fake_records)
        self.assertEqual(fake_records, procesed)

    def test_process_wgs_master_scaffolds_contig(self):
        """Test process_wgs_master_scaffolds() gets contigs from NCBI if there is a CONTIG annotation"""
        mock('utils.fix_supercontig_record', tracker=self.tt, returns_func=lambda x: [x])
        fake_records = [
            FakeRecord(annotations={'contig': 'foo'}),
            FakeRecord(annotations={'contig': 'bar'})
        ]
        fake_records[0].seq = UnknownSeq(23)
        fake_records[1].seq = UnknownSeq(42)
        trace = """Called utils.fix_supercontig_record(
    <antismash.test.test_utils.FakeRecord object at ...>)
Called utils.fix_supercontig_record(
    <antismash.test.test_utils.FakeRecord object at ...>)
"""
        procesed = utils.process_wgs_master_scaffolds(fake_records)
        self.assertEqual(fake_records, procesed)
        assert_same_trace(self.tt, trace)

    def test_process_wgs_master_scaffolds_wgs(self):
        """Test process_wgs_master_scaffolds() gets contigs from NCBI if there are WGS/WGS_SCAFLD annotations"""
        mock('utils.fix_wgs_master_record', tracker=self.tt, returns_func=lambda x: [x])
        fake_records = [
            FakeRecord(annotations={'wgs': 'foo'}),
            FakeRecord(annotations={'wgs_scafld': 'bar'})
        ]
        fake_records[0].seq = UnknownSeq(23)
        fake_records[1].seq = UnknownSeq(42)
        trace = """Called utils.fix_wgs_master_record(
    <antismash.test.test_utils.FakeRecord object at ...>)
Called utils.fix_wgs_master_record(
    <antismash.test.test_utils.FakeRecord object at ...>)
"""
        procesed = utils.process_wgs_master_scaffolds(fake_records)
        self.assertEqual(fake_records, procesed)
        assert_same_trace(self.tt, trace)

    def test_process_wgs_master_scaffolds_contig_with_seq(self):
        """Test process_wgs_master_scaffolds() doesn't get contigs if there is a CONTIG ann. but also a sequence"""
        fake_records = [
            FakeRecord(annotations={'contig': 'foo'}),
            FakeRecord(annotations={'contig': 'bar'})
        ]
        procesed = utils.process_wgs_master_scaffolds(fake_records)
        self.assertEqual(fake_records, procesed)

    def test_process_wgs_master_scaffolds_wgs_with_seq(self):
        """Test process_wgs_master_scaffolds() doesn't get contigs if there are WGS/WGS_SCAFLDs but also sequences"""
        mock('utils.fix_wgs_master_record', tracker=self.tt, returns_func=lambda x: [x])
        fake_records = [
            FakeRecord(annotations={'wgs': 'foo'}),
            FakeRecord(annotations={'wgs_scafld': 'bar'})
        ]

        procesed = utils.process_wgs_master_scaffolds(fake_records)
        self.assertEqual(fake_records, procesed)
