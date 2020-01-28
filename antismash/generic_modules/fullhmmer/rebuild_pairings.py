#!/usr/bin/env python
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

from __future__ import print_function

"""
(Re)builds an importable python dictionary from a PFAM database, mapping PFAM
name to accession (e.g. CbiA -> PF01656).

By default looks for database in the same directory named 'Pfam-A.hmm' and
outputs to a python file 'name2pfamid.py'.

"""

HEADER = """# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

name_to_pfamid = {
"""

FOOTER = """
}
"""

def fetch_value(line):
    value = line.split()[1].strip()
    assert value, "empty value in line: %s" % line
    return value

def rebuild_pairings(database_filename="Pfam-A.hmm", output_filename="name2pfamid.py"):
    names = []
    accessions = []
    with open(database_filename, "r") as database:
        for line in database:
            if line.startswith("NAME"):
                names.append(fetch_value(line))
            elif line.startswith("ACC"):
                # strip any trailing info from accession: e.g. PF01656.18 -> PF01656
                accessions.append(fetch_value(line).split(".", 1)[0])
    assert len(names) == len(accessions), "some pairings incomplete"
    lines = []
    for name, acc in zip(names, accessions):
        lines.append('    "{}" : "{}"'.format(name, acc))
    with open(output_filename, "w") as pairings:
        pairings.write(HEADER)
        pairings.write(",\n".join(lines))
        pairings.write(FOOTER)
    return len(lines)

if __name__ == "__main__":
    count = rebuild_pairings()
    print("rebuilt with %d pairings" % count)
