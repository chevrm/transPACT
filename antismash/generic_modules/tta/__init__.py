# vim: set fileencoding=utf-8 :
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

"""Identify TTA codons in BGCs"""

import logging
from Bio.SeqFeature import SeqFeature, FeatureLocation
from antismash import utils


def check_prereqs(options):
    """Check for prerequesites"""
    # No external dependencies
    return []

def detect(seq_record, options):
    """Detect TTA codons"""
    logging.info("Detecting TTA codons")
    cds_features = utils.get_withincluster_cds_features(seq_record)
    for feature in cds_features:
        sequence = feature.extract(seq_record.seq)
        for i in xrange(0, len(sequence), 3):
            codon = sequence[i:i+3].lower()
            if "tta" == codon:
                seq_record.features.append(_create_tta_feature(feature, i))


def _create_tta_feature(feature, offset):
    """Create a misc_feature entry for a TTA codon on a given feature"""
    if feature.strand == 1:
        start = feature.location.start + offset
        end = start + 3
    else:
        end = feature.location.end - offset
        start = end - 3

    loc = FeatureLocation(start, end, feature.strand)

    qualifiers = {
        "note": ["tta leucine codon, possible target for bldA regulation"],
        "tool": ["antiSMASH"]
    }

    tta_feature = SeqFeature(loc, type="misc_feature", qualifiers=qualifiers)
    return tta_feature
