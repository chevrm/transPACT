# vim: set fileencoding=utf-8 :
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
import logging
from os import path
from antismash import utils


def resolve_predicat_domain_specificity(predicat_domain):
    id_line = find_full_id_line(predicat_domain)
    return get_specificity_from_id_line(id_line)


def find_full_id_line(predicat_domain):
    from antismash.specific_modules.nrpspks import predicat_id_lines

    for line in predicat_id_lines:
        if line.startswith(predicat_domain):
            return line
        else:
            mangled_line = line.replace('(', '-').replace(')', '-')
            if mangled_line.startswith(predicat_domain):
                return line
    logging.debug("No match for %r", predicat_domain)
    return "no_match_n/a"


def get_specificity_from_id_line(id_line):
    parts = id_line.split('_')
    return parts[-1]


def load_id_lines():
    sandpuma_dir = utils.get_full_path(__file__, 'sandpuma')
    fasta_file = path.join(sandpuma_dir, 'flat', 'fullset0_smiles.faa')

    id_lines = []

    with open(fasta_file, 'r') as fh:
        for line in fh:
            if not line.startswith(">"):
                continue

            id_lines.append(line.strip().lstrip(">"))

    return id_lines
