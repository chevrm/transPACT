"Test suite for hmm_detection module"

try:
    import unittest2
except ImportError:
    import unittest as unittest2
from minimock import mock, restore, TraceTracker, assert_same_trace
from argparse import Namespace

from antismash import utils
from antismash import config
from antismash.generic_modules import hmm_detection
from Bio.SeqFeature import SeqFeature, FeatureLocation

class FakeRecord(object):
    "class for generating a seq_record like data structure"
    def __init__(self, features=None, seq='FAKESEQ'):
        if features is None:
            features = []
        self.features = features
        self.seq = FakeSeq(seq)
    def __len__(self):
        return 1000

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

    def translate(self, to_stop):
        return self.seq

    def __str__(self):
        return self.seq


class FakeHSP(object):
    "class for generating a HSP like datastructure"
    def __init__(self, query_id, hit_id, hit_start, hit_end, bitscore, evalue):
        self.query_id = query_id
        self.hit_id = hit_id
        self.hit_start = hit_start
        self.hit_end = hit_end
        self.bitscore = bitscore
        self.evalue = evalue


class HmmDetectionTest(unittest2.TestCase):
    def setUp(self):
        self.config = Namespace()
        config.set_config(self.config)
        self.results_by_id = {
            "GENE_1" : [
                FakeHSP("modelA", "GENE_1", 0, 10, 50, 0),
                FakeHSP("modelB", "GENE_1", 0, 10, 50, 0)
            ],
            "GENE_2" : [
                FakeHSP("modelA", "GENE_1", 0, 10, 50, 0),
                FakeHSP("modelB", "GENE_1", 0, 10, 50, 0)
            ],
            "GENE_3" : [
                FakeHSP("modelA", "GENE_1", 0, 10, 50, 0),
                FakeHSP("modelB", "GENE_1", 0, 10, 50, 0)
            ],
            "GENE_4" : [
                FakeHSP("modelA", "GENE_1", 0, 10, 50, 0),
                FakeHSP("modelB", "GENE_1", 0, 10, 50, 0)
            ],
            "GENE_5" : [
                FakeHSP("modelA", "GENE_1", 0, 10, 50, 0),
                FakeHSP("modelB", "GENE_1", 0, 10, 50, 0)
            ]
        }
        self.feature_by_id = {
            "GENE_1" : FakeFeature("CDS", FeatureLocation(0, 30), {"locus_tag": ["GENE_1"]}),
            "GENE_2" : FakeFeature("CDS", FeatureLocation(30, 50), {"locus_tag": ["GENE_2"]}),
            "GENE_3" : FakeFeature("CDS", FeatureLocation(70, 90), {"locus_tag": ["GENE_3"]}),
            "GENE_X" : FakeFeature("CDS", FeatureLocation(95, 100), {"locus_tag": ["GENE_X"]}),
            "GENE_4" : FakeFeature("CDS", FeatureLocation(120, 130), {"locus_tag": ["GENE_4"]}),
            "GENE_5" : FakeFeature("CDS", FeatureLocation(130, 150), {"locus_tag": ["GENE_5"]})
        }
        self.rulesdict = {
            "MetaboliteA" : ("modelA", 10, 5),
            "MetaboliteB" : ("(modelA & modelB)", 10, 5),
            "MetaboliteC" : ("cluster(modelA,modelB)", 10, 5),
            "MetaboliteD" : ("minimum(3,[modelA,modelB], [modelA])", 20, 5),
            "Metabolite0" : ("modelC", 1, 3),
            "Metabolite1" : ("modelC", 1, 3)
        }
        self.features = []
        for gene_id in self.feature_by_id:
            self.features.append(self.feature_by_id[gene_id])
        self.record = FakeRecord(self.features)

    def test_apply_cluster_rules(self):
        enabled_clustertypes = list(set(self.rulesdict.keys()))
        detected_types = hmm_detection.apply_cluster_rules(self.results_by_id, self.feature_by_id, enabled_clustertypes, self.rulesdict, utils.get_overlaps_table(self.record))
        for gid in detected_types:
            detected_types[gid] = set(detected_types[gid].split("-"))
        expected_types = {
            "GENE_1" : set(["MetaboliteA", "MetaboliteB", "MetaboliteC", "MetaboliteD"]),
            "GENE_2" : set(["MetaboliteA", "MetaboliteB", "MetaboliteC", "MetaboliteD"]),
            "GENE_3" : set(["MetaboliteA", "MetaboliteB", "MetaboliteC", "MetaboliteD"]),
            "GENE_4" : set(["MetaboliteA", "MetaboliteB", "MetaboliteC"]),
            "GENE_5" : set(["MetaboliteA", "MetaboliteB", "MetaboliteC"])
        }
        self.assertEqual(detected_types, expected_types, msg = "\nResult : %s\nExpected : %s" % (detected_types, expected_types))

    def test_find_clusters(self):
        i = 0
        nseqdict = {"Metabolite0" : "?", "Metabolite1" : "?"}
        self.config.next_clusternr = 1
        for gene_id in self.feature_by_id:
            if gene_id != "GENE_X":
                clustertype = "Metabolite%d" % (i % 2)
                hmm_detection._update_sec_met_entry(self.feature_by_id[gene_id], self.results_by_id[gene_id], clustertype, nseqdict)
                i += 1
        hmm_detection.find_clusters(self.record, self.rulesdict)
        result_clusters = [sorted([utils.get_gene_id(f) for f in utils.get_cluster_cds_features(feature,self.record)]) for feature in utils.get_cluster_features(self.record)]
        expected_clusters = [
            ["GENE_1", "GENE_2"],
            ["GENE_3"],
            ["GENE_4", "GENE_5"]
        ]
        self.assertEqual(result_clusters, expected_clusters, msg = "\nResult : %s\nExpected : %s" % (result_clusters, expected_clusters))
