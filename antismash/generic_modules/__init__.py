class Signature(object):
    """Secondary metabolite signature"""
    def __init__(self, name, _type, description, cutoff, path):
        self.name = name
        self.type = _type
        self.description = description
        self.cutoff = cutoff
        self.path = path

from antismash.generic_modules import (
        active_site_finder,
        fullhmmer,
#        fullhmmer_dblookup,
        genefinding,
        hmm_detection,
        clusterblast,
        subclusterblast,
        knownclusterblast,
        smcogs,
        cassis,
        tta,
        gff_parser,
    )


def check_prereqs(options):
    failure_msgs = []

    if options.full_hmmer or options.inclusive:
        failure_msgs.extend(fullhmmer.check_prereqs(options))

    failure_msgs.extend(genefinding.check_prereqs(options))
    failure_msgs.extend(hmm_detection.check_prereqs())

    if options.run_asf:
        failure_msgs.extend(active_site_finder.check_prereqs(options))

    if options.tta:
        failure_msgs.extend(tta.check_prereqs(options))

    if options.smcogs:
        failure_msgs.extend(smcogs.check_prereqs(options))

    if options.clusterblast:
        failure_msgs.extend(clusterblast.check_prereqs(options))

    if options.subclusterblast:
        failure_msgs.extend(subclusterblast.check_prereqs(options))

    if options.knownclusterblast:
        failure_msgs.extend(knownclusterblast.check_prereqs(options))

    if options.cassis:
        failure_msgs.extend(cassis.check_prereqs(options))

    return failure_msgs
