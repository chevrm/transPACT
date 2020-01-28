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

import logging
from antismash import utils
import os
from os import path
import shutil
from helperlibs.wrappers.io import TemporaryDirectory

def extract_nrps_genes(pksnrpscoregenes, domaindict, seq_record, extra_aa=0):
    nrpsnames = [] 
    nrpsseqs = []
    for feature in pksnrpscoregenes:
        locus = utils.get_gene_id(feature)
        domaindetails = domaindict[locus]
        nr = 0
        for tab in domaindetails:
            if tab[0] == "AMP-binding" or tab[0] == "A-OX":
                nr += 1
                start = int(tab[1])
                end = int(tab[2]) + extra_aa
                seq = str(utils.get_aa_sequence(feature))[start:end]
                name = locus + "_A" + str(nr)
                nrpsnames.append(name)
                nrpsseqs.append(seq)
    return nrpsnames, nrpsseqs

def run_sandpuma(seq_record, nrpsnames, nrpsseqs, options):
    """Run SANDPUMA on the set of NRPS sequences from this genome"""

    nrpspredictor_output = "ctg" + str(options.record_idx) + "_nrpspredictor3_svm.txt"
    individual_predictions = "ctg" + str(options.record_idx) + "_ind.res.tsv"
    percentage_identities = "ctg" + str(options.record_idx) + "_pid.res.tsv"
    sandpuma_predictions = "ctg" + str(options.record_idx) + "_sandpuma.tsv"
    ensemble_predictions = "ctg" + str(options.record_idx) + "_ens.res.tsv"

    # In debug mode, simply copy over previous predictions.
    if options.dbgsandpuma != '':
        shutil.copy(path.join(options.dbgsandpuma, nrpspredictor_output),
                    path.join(options.raw_predictions_outputfolder, nrpspredictor_output))
        shutil.copy(path.join(options.dbgsandpuma, individual_predictions),
                    path.join(options.raw_predictions_outputfolder, individual_predictions))
        shutil.copy(path.join(options.dbgsandpuma, percentage_identities),
                    path.join(options.raw_predictions_outputfolder, percentage_identities))
        shutil.copy(path.join(options.dbgsandpuma, sandpuma_predictions),
                    path.join(options.raw_predictions_outputfolder, sandpuma_predictions))
        shutil.copy(path.join(options.dbgsandpuma, ensemble_predictions),
                    path.join(options.raw_predictions_outputfolder, ensemble_predictions))
        return

    
    #NRPSPredictor: extract AMP-binding + 120 residues N-terminal of this domain, extract 8 Angstrom residues and insert this into NRPSPredictor
    sandpumadir = utils.get_full_path(__file__, "sandpuma")
    with TemporaryDirectory(change=True):
        #Extract A domains from the NRPS sequences and write to FASTA file
        nrpsseqs_file = "input_adomains.fasta"
        utils.writefasta(nrpsnames, nrpsseqs, nrpsseqs_file)
        #Run SANDPUMA on the FASTA file
        sandpuma_command = [sandpumadir + os.sep + 'predictnrps_nodep_par.sh', 'input_adomains.fasta', sandpumadir, str(options.cpus)]
        err = utils.execute(sandpuma_command)[1]
        if err != '':
            logging.error('Running SANDPUMA gave an error')
            raise RuntimeError("Sandpuma failed to run: %s" % err)
        #Copy SANDPUMA (including NRPSPredictor2) results and move back to original directory
        shutil.move("query.rep", options.raw_predictions_outputfolder + os.sep + nrpspredictor_output)
        shutil.move("ind.res.tsv", options.raw_predictions_outputfolder + os.sep + individual_predictions)
        shutil.move("pid.res.tsv", options.raw_predictions_outputfolder + os.sep + percentage_identities)
        shutil.move("sandpuma.tsv", options.raw_predictions_outputfolder + os.sep + sandpuma_predictions)
        shutil.move("ens.res.tsv", options.raw_predictions_outputfolder + os.sep + ensemble_predictions)
