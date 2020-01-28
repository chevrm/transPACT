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

import os
from os import path
from collections import defaultdict
from antismash import utils


def parse_nrps_preds(options, pksnrpsvars):
    pksnrpsvars.sandpuma_allmethods = {}
    pksnrpsvars.sandpuma_pid = {}
    pksnrpsvars.sandpuma_snn = {}
    pksnrpsvars.sandpuma_res = defaultdict(lambda: "N/A")

    individual_predictions = "ctg" + str(options.record_idx) + "_ind.res.tsv"
    percentage_identities = "ctg" + str(options.record_idx) + "_pid.res.tsv"
    sandpuma_predictions = "ctg" + str(options.record_idx) + "_sandpuma.tsv"

    if len(pksnrpsvars.nrpsnames) > 0:
        with open(path.join(options.raw_predictions_outputfolder, individual_predictions), 'r') as fh:
            for line in fh:
                line = line.strip()
                if line == '':
                    break

                tabs = line.split("\t")
                if tabs[0] not in pksnrpsvars.sandpuma_allmethods.keys():
                    pksnrpsvars.sandpuma_allmethods[tabs[0]] = {}
                pksnrpsvars.sandpuma_allmethods[tabs[0]][tabs[1]] = tabs[2]
                if len(tabs) > 3:
                    pksnrpsvars.sandpuma_snn[tabs[0]] = tabs[3]

        with open(path.join(options.raw_predictions_outputfolder, percentage_identities), 'r') as fh:
            for line in fh:
                line = line.strip()
                if line == '':
                    break

                tabs = line.split("\t")
                pksnrpsvars.sandpuma_pid[tabs[0]] = tabs[1]

        with open(path.join(options.raw_predictions_outputfolder, sandpuma_predictions), 'r') as fh:
            for line in fh:
                line = line.strip()
                if line == '':
                    break

                tabs = line.split("\t")
                pksnrpsvars.sandpuma_res[tabs[0]] = tabs[2]


def parse_pks_preds(options, pksnrpsvars):
    pksnrpsvars.minowa_pks_preds_details = {}
    pksnrpsvars.minowa_pks_preds = {}
    pksnrpsvars.pks_code_preds = {}
    pksnrpsvars.pks_code_preds_details = {}
    pksnrpsvars.minowa_cal_preds = {}
    pksnrpsvars.minowa_cal_preds_details = {}
    substratetransdict = {'Malonyl-CoA': 'mal', 'Methylmalonyl-CoA': 'mmal', 'Methoxymalonyl-CoA': 'mxmal',
                          'Ethylmalonyl-CoA': 'emal', 'Isobutyryl-CoA': 'isobut', '2-Methylbutyryl-CoA': '2metbut',
                          'trans-1,2-CPDA': 'trans-1,2-CPDA', 'Acetyl-CoA': 'Acetyl-CoA', 'Benzoyl-CoA': 'benz',
                          'Propionyl-CoA': 'prop', '3-Methylbutyryl-CoA': '3metbut', 'Ethylmalonyl-CoA': 'Ethyl_mal',
                          'CE-Malonyl-CoA': 'cemal', '2-Rhyd-Malonyl-CoA': '2Rhydmal', 'CHC-CoA': 'CHC-CoA',
                          'inactive': 'inactive'}
    if len(pksnrpsvars.pksnames) > 0:
        minowa_at_file = open(options.raw_predictions_outputfolder + os.sep + "ctg" + str(
            options.record_idx) + "_minowa_pkspredoutput.txt", "r")
        minowa_at_file = minowa_at_file.read()
        minowa_at_file = minowa_at_file.replace("\r", "\n")
        parts = minowa_at_file.split("\\\\\n")[1:]
        for i in parts:
            partlines = i.split("\n")
            acc = partlines[0]
            if substratetransdict.has_key(partlines[2].split("\t")[0]):
                tophit = substratetransdict[partlines[2].split("\t")[0]]
            else:
                tophit = "pk"
            pksnrpsvars.minowa_pks_preds[acc] = tophit
            pksnrpsvars.minowa_pks_preds_details[
                acc] = "<b>Minowa HMM method AT-domain<br>Substrate specificity prediction top hits:</b><br>\n" + \
                       partlines[1] + "<br>\n" + partlines[2] + "<br>\n" + partlines[3] + "<br>\n" + partlines[
                           4] + "<br><br>\n\n"
        pkssignaturefile = open(
            options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_pkssignatures.txt", "r")
        pkssignaturefile = pkssignaturefile.read()
        pkssignaturefile = pkssignaturefile.replace("\r", "\n")
        parts = pkssignaturefile.split("//\n")[1:]
        for i in parts:
            partlines = i.split("\n")
            partlines2 = []
            for j in partlines:
                if j != "":
                    partlines2.append(j)
            partlines = partlines2
            acc = partlines[0].split("\t")[0]
            if len(partlines) > 2:
                tophit = (partlines[1].split("\t")[0]).split("__")[1]
                pksnrpsvars.pks_code_preds[acc] = tophit
                codes = []
                prots = []
                scores = []
                for i in partlines[1:4]:
                    codes.append(i.split("\t")[0])
                    prot = i.split("\t")[1]
                    prot = prot.replace("_AT", " (AT")
                    prot = prot.replace("__", "): ")
                    prots.append(prot)
                    scores.append(i.split("\t")[2])
                if len(prots) >= 3:
                    pksnrpsvars.pks_code_preds_details[
                        acc] = "<b>PKS Active Site Signature method<br>AT-domain substrate specificity prediction top hits:</b><br>\nCode:" + \
                               partlines[0].split("\t")[1] + "<br>\n" + codes[0] + " - " + prots[0] + " : (" + scores[
                                   0] + "% identity)<br>\n" + codes[1] + " - " + prots[1] + " : (" + scores[
                                   1] + "% identity)<br>\n" + codes[2] + " - " + prots[2] + " : (" + scores[
                                   2] + "% identity)<br><br>\n\n"
                elif len(prots) == 2:
                    pksnrpsvars.pks_code_preds_details[
                        acc] = "<b>PKS Active Site Signature method<br>AT-domain substrate specificity prediction top hits:</b><br>\nCode:" + \
                               partlines[0].split("\t")[1] + "<br>\n" + codes[0] + " - " + prots[0] + " : (" + scores[
                                   0] + "% identity)<br>\n" + codes[1] + " - " + prots[1] + " : (" + scores[
                                   1] + "% identity)<br><br>\n\n"
                elif len(prots) == 1:
                    pksnrpsvars.pks_code_preds_details[
                        acc] = "<b>PKS Active Site Signature method<br>AT-domain substrate specificity prediction top hits:</b><br>\nCode:" + \
                               partlines[0].split("\t")[1] + "<br>\n" + codes[0] + " - " + prots[0] + " : (" + scores[
                                   0] + "% identity)<br><br>\n\n"
            else:
                pksnrpsvars.pks_code_preds[acc] = "N/A"
                pksnrpsvars.pks_code_preds_details[
                    acc] = "<b>PKS Active Site Signature method<br>No AT-domain substrate specificity prediction hits above 40% identity.<br>\n\n"
    if len(pksnrpsvars.calnames) > 0:
        minowa_cal_file = open(options.raw_predictions_outputfolder + os.sep + "ctg" + str(
            options.record_idx) + "_minowa_calpredoutput.txt", "r")
        minowa_cal_file = minowa_cal_file.read()
        minowa_cal_file = minowa_cal_file.replace("\r", "\n")
        parts = minowa_cal_file.split("\\\\\n")[1:]
        for i in parts:
            partlines = i.split("\n")
            acc = partlines[0]
            tophit = partlines[2].split("\t")[0]
            pksnrpsvars.minowa_cal_preds[acc] = tophit
            pksnrpsvars.minowa_cal_preds_details[
                acc] = "<b>Minowa HMM method<br>CAL-domain substrate specificity prediction top hits:</b><br>\n" + \
                       partlines[1] + "<br>\n" + partlines[2] + "<br>\n" + partlines[3] + "<br>\n" + partlines[
                           4] + "<br><br>\n\n"


def parse_kr_activity_preds(options, pksnrpsvars):
    pksnrpsvars.kr_activity_preds = {}
    pksnrpsvars.kr_stereo_preds = {}
    if len(pksnrpsvars.krnames) > 0:
        krfile = open(
            options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_krpredoutput.txt", "r")
        krfile = krfile.read()
        krfile = krfile.replace("\r", "\n")
        krlines = krfile.split("\n")[:-1]
        for i in krlines:
            tabs = i.split("\t")
            pksnrpsvars.kr_activity_preds[tabs[0]] = tabs[1]
            pksnrpsvars.kr_stereo_preds[tabs[0]] = tabs[2]


def calculate_consensus_prediction(pksnrpsvars, seq_record):
    # Combine substrate specificity predictions into consensus prediction
    pksnrpsvars.consensuspreds = {}
    pksnrpsvars.consensuspreds_transat = {}
    available_smiles_parts = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'PHE', 'TRP', 'SER', 'THR', 'ASN', 'GLN',
                              'TYR', 'CYS', 'LYS', 'ARG',
                              'HIS', 'ASP', 'GLU', 'MPRO', 'ORN', 'PGLY', 'DAB', 'BALA', 'AEO', 'DHA', 'PIP', 'BMT',
                              'gly', 'ala', 'val', 'leu', 'ile', 'met', 'pro', 'phe', 'trp', 'ser',
                              'thr', 'asn', 'gln', 'tyr', 'cys', 'lys', 'arg', 'his', 'asp', 'glu', 'aaa', 'mpro',
                              'dhb', '2hiva', 'orn', 'pgly', 'dab', 'bala', 'aeo', '4mha', 'pico', 'phg',
                              'dha', 'scy', 'pip', 'bmt', 'adds', 'aad', 'abu', 'hiv', 'dhpg', 'bht', '3-me-glu',
                              '4pPro', 'ala-b', 'ala-d', 'dht', 'Sal', 'tcl', 'lys-b', 'hpg', 'hyv-d',
                              'iva', 'vol', 'mal', 'mmal', 'ohmal', 'redmal', 'mxmal', 'emal', 'nrp', 'pk', 'Gly',
                              'Ala', 'Val', 'Leu', 'Ile', 'Met', 'Pro', 'Phe', 'Trp', 'Ser', 'Thr', 'Asn', 'Gln', 'Tyr',
                              'Cys', 'Lys', 'Arg', 'His', 'Asp', 'Glu', 'Mpro', '23Dhb', '34Dhb', '2Hiva', 'Orn',
                              'Pgly', 'Dab', 'Bala', 'Aeo', '4Mha', 'Pico', 'Aaa', 'Dha', 'Scy', 'Pip',
                              'Bmt', 'Adds', 'DHpg', 'DHB', 'nrp', 'pk']

    # Extracting gene cluster type (e.g., "transatpks")
    for f in utils.get_cluster_features(seq_record):
        cluster_info = f.qualifiers

    for feature in pksnrpsvars.pksnrpscoregenes:
        locus = utils.get_gene_id(feature)
        nra = 0
        nrat = 0
        nrcal = 0
        nrtransat = 0
        j = pksnrpsvars.domaindict[locus]

        for k in j:
            if 'transatpks' not in cluster_info['product'][0]:
                if k[0] == "PKS_AT":
                    nrat += 1
                    preds = []
                    preds.append(pksnrpsvars.minowa_pks_preds[locus + "_AT" + str(nrat)])
                    preds.append(pksnrpsvars.pks_code_preds[locus + "_AT" + str(nrat)])
                    cpred = "n"
                    for l in preds:
                        if preds.count(l) > 1:
                            if l in available_smiles_parts:
                                pksnrpsvars.consensuspreds[locus + "_AT" + str(nrat)] = l
                            else:
                                pksnrpsvars.consensuspreds[locus + "_AT" + str(nrat)] = "pk"
                            cpred = "y"
                    if cpred == "n":
                        pksnrpsvars.consensuspreds[locus + "_AT" + str(nrat)] = "pk"
            elif 'transatpks' in cluster_info['product'][0]:
                if k[0] == "PKS_AT":
                    nrat += 1
                    preds = []
                    preds.append(pksnrpsvars.minowa_pks_preds[locus + "_AT" + str(nrat)])
                    preds.append(pksnrpsvars.pks_code_preds[locus + "_AT" + str(nrat)])
                    cpred = "n"

                    # Only for the writing purpose in sec_record (i.e., trans-AT)
                    for l in preds:
                        if preds.count(l) > 1:
                            if l in available_smiles_parts:
                                pksnrpsvars.consensuspreds_transat[locus + "_AT" + str(nrat)] = l
                            else:
                                pksnrpsvars.consensuspreds_transat[locus + "_AT" + str(nrat)] = "pk"
                            cpred = "y"
                    if cpred == "n":
                        pksnrpsvars.consensuspreds_transat[locus + "_AT" + str(nrat)] = "pk"
                # For chemical display purpose for chemicals from trans-AT PKS gene cluster
                # mal is always assumed for trans-AT
                if k[0] == "PKS_KS":
                    nrtransat += 1
                    pksnrpsvars.consensuspreds[locus + "_KS" + str(nrtransat)] = "mal"
                    cpred = "y"
            if k[0] == "AMP-binding" or k[0] == "A-OX":
                nra += 1
                if pksnrpsvars.sandpuma_res[locus + "_A" + str(nra)] == "no_call":
                    pksnrpsvars.consensuspreds[locus + "_A" + str(nra)] = "nrp"
                else:
                    pksnrpsvars.consensuspreds[locus + "_A" + str(nra)] = pksnrpsvars.sandpuma_res[
                        locus + "_A" + str(nra)]
            if k[0] == "CAL_domain":
                nrcal += 1
                if pksnrpsvars.minowa_cal_preds[locus + "_CAL" + str(nrcal)] in available_smiles_parts:
                    pksnrpsvars.consensuspreds[locus + "_CAL" + str(nrcal)] = pksnrpsvars.minowa_cal_preds[
                        locus + "_CAL" + str(nrcal)]
                else:
                    pksnrpsvars.consensuspreds[locus + "_CAL" + str(nrcal)] = "pk"


def generate_domainnamesdict(pksnrpsvars):
    pksnrpsvars.domainnamesdict = {}
    for feature in pksnrpsvars.pksnrpscoregenes:
        locus = utils.get_gene_id(feature)
        j = pksnrpsvars.domaindict[locus]
        domainnames = []
        for k in j:
            domainnames.append(k[0])
        pksnrpsvars.domainnamesdict[locus] = domainnames


def find_duplicate_position(domainList, item):
    start_at = -1
    locs = []
    while True:
        try:
            loc = domainList.index(item, start_at + 1)
        except ValueError:
            break
        else:
            locs.append(loc)
            start_at = loc
    return locs

def update_prediction(locus, preds, target, target_list, lists, mappings):
    assert len(lists) == len(mappings)
    for idx in range(len(target_list)):
        key = locus + target + str(idx + 1)
        for sublist, mapping in zip(lists, mappings):
            for position in sublist:
                if not target_list[idx] < position:
                    continue
                if idx + 1 <= len(target_list) - 1 and not position < target_list[idx + 1]:
                    continue
                current = preds[key]
                preds[key] = mapping.get(current, current)

def insert_modified_monomers(pksnrpsvars, seq_record):

    # Extracting gene cluster type (e.g., "transatpks")
    for f in utils.get_cluster_features(seq_record):
        cluster_info = f.qualifiers

    # pksnrpsvars.domainnamesdict = {'CRYAR_RS43165': ['PKS_KS', 'PKS_AT',...]}
    # Get a unique set of genes having ATs
    locusTag_domain = sorted(set(pksnrpsvars.domainnamesdict.keys()))
    # for modifications, e.g. mal -> ohmal
    # must be the same length and ordering as the lists generation below
    mappings = [{"mal":"ohmal", "mmal":"ohmmal", "mxmal":"ohmxmal", "emal":"ohemal"}, # KR
                {"ohmal":"ccmal", "ohmmal":"ccmmal", "ohmxmal":"ccmxmal", "ohemal":"ccemal"}, # DH
                {"ccmal":"redmal", "ccmmal":"redmmal", "ccmxmal":"redmxmal", "ccemal":"redemal"}] # ER

    for locusTag in locusTag_domain:
        at_list = find_duplicate_position(pksnrpsvars.domainnamesdict[locusTag], 'PKS_AT')
        # For transatpks
        ks_list = find_duplicate_position(pksnrpsvars.domainnamesdict[locusTag], 'PKS_KS')

        lists = [find_duplicate_position(pksnrpsvars.domainnamesdict[locusTag], 'PKS_KR'),
                 find_duplicate_position(pksnrpsvars.domainnamesdict[locusTag], 'PKS_DH'),
                 find_duplicate_position(pksnrpsvars.domainnamesdict[locusTag], 'PKS_ER')]

        if 'transatpks' not in cluster_info['product'][0]:
            update_prediction(locusTag, pksnrpsvars.consensuspreds, "_AT", at_list, lists, mappings)
        else:
            update_prediction(locusTag, pksnrpsvars.consensuspreds, "_KS", ks_list, lists, mappings)


def write_substr_spec_preds_to_html(options, pksnrpsvars):
    # Make output folder for storing HTML files with substrate specificities
    options.substrspecsfolder = options.full_outputfolder_path + os.sep + "html"
    if not os.path.exists(options.substrspecsfolder):
        os.mkdir(options.substrspecsfolder)
    # Write all prediction details to HTML files for each gene to be used as pop-up window
    for feature in pksnrpsvars.pksnrpscoregenes:
        locus = utils.get_gene_id(feature)
        if "PKS_AT" in pksnrpsvars.domainnamesdict[locus] or "AMP-binding" in pksnrpsvars.domainnamesdict[
            locus] or "A-OX" in pksnrpsvars.domainnamesdict[locus] or "CAL_domain" in pksnrpsvars.domainnamesdict[
            locus]:
            j = pksnrpsvars.domaindict[locus]
            nrat = 0
            nra = 0
            nrcal = 0
            for k in j:
                if k[0] == "PKS_AT":
                    nrat += 1
                    domainname = locus + "_AT" + str(nrat)
                    htmloutfile = open(options.substrspecsfolder + os.sep + domainname + ".html", "w")
                    htmloutfile.write(
                        '<html>\n<head>\n<title>Prediction details</title>\n<STYLE type="text/css">\nbody{\n  text-align:left;\n  background-color:white;\n  font-family: Tahoma, sans-serif;\n  font-size: 0.8em;\n  color: #810E15;\n}\n</STYLE>\n</head>\n<body>')
                    htmloutfile.write(pksnrpsvars.minowa_pks_preds_details[domainname])
                    htmloutfile.write(pksnrpsvars.pks_code_preds_details[domainname])
                    htmloutfile.write(
                        "<b><i>Consensus Predictions: " + pksnrpsvars.consensuspreds[domainname] + "</b></i>")
                    htmloutfile.write('\n</body>\n</html>')
                    htmloutfile.close()
                if k[0] == "AMP-binding" or k[0] == "A-OX":
                    nra += 1
                    domainname = locus + "_A" + str(nra)
                    htmloutfile = open(options.substrspecsfolder + os.sep + domainname + ".html", "w")
                    htmloutfile.write(
                        '<html>\n<head>\n<title>Prediction details</title>\n<STYLE type="text/css">\nbody{\n  text-align:left;\n  background-color:white;\n  font-family: Tahoma, sans-serif;\n  font-size: 0.8em;\n  color: #810E15;\n}\n</STYLE>\n</head>\n<body>')
                    htmloutfile.write(pksnrpsvars.nrps_svm_preds_details[domainname])
                    htmloutfile.write(pksnrpsvars.nrps_code_preds_details[domainname])
                    htmloutfile.write(pksnrpsvars.minowa_nrps_preds_details[domainname])
                    htmloutfile.write(
                        "<b><i>Consensus Prediction: '" + pksnrpsvars.consensuspreds[domainname] + "'</b></i>")
                    htmloutfile.write('\n</body>\n</html>')
                    htmloutfile.close()
                if k[0] == "CAL_domain":
                    nrcal += 1
                    domainname = locus + "_CAL" + str(nrcal)
                    htmloutfile = open(options.substrspecsfolder + os.sep + domainname + ".html", "w")
                    htmloutfile.write(
                        '<html>\n<head>\n<title>Prediction details</title>\n<STYLE type="text/css">\nbody{\n  text-align:left;\n  background-color:white;\n  font-family: Tahoma, sans-serif;\n  font-size: 0.8em;\n  color: #810E15;\n}\n</STYLE>\n</head>\n<body>')
                    htmloutfile.write(pksnrpsvars.minowa_cal_preds_details[domainname])
                    htmloutfile.write('\n</body>\n</html>')
                    htmloutfile.close()
