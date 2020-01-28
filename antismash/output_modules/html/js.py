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

from antismash import utils
import logging
import re
import os
from glob import glob
from antismash.output_modules.html.generate_html_table import generate_html_table

def convert_records(seq_records, options):
    records = []
    annotations = load_cog_annotations()
    for srec in seq_records:
        records.append(convert_record(srec, annotations, options))
    return records

def convert_record(record, annotations, options):
    """Convert a SeqRecord to JSON"""
    js_rec = {}
    js_rec['seq_id'] = record.id
    if "extrarecord" in options:
        if options.extrarecord.has_key(record.id):
            if "extradata" in options.extrarecord[record.id]:
                if options.extrarecord[record.id].extradata.has_key("orig_id"):
                    js_rec['orig_id'] = options.extrarecord[record.id].extradata["orig_id"]
    if not js_rec.has_key('orig_id'):
        js_rec['orig_id'] = ""
    js_rec['clusters'] = convert_clusters(record, annotations, options)

    return js_rec

def convert_clusters(record, annotations, options):
    """Convert cluster SeqFeatures to JSON"""
    js_clusters = []
    for cluster in utils.get_cluster_features(record):
        features = utils.get_cluster_cds_features(cluster, record)
        borders = utils.get_cluster_cluster_border_features(cluster, record)

        tta_codons = []
        all_misc_features = utils.get_all_features_of_type(record, 'misc_feature')
        for feature in all_misc_features:
            if not utils.features_overlap(cluster, feature):
                continue
            if 'note' not in feature.qualifiers:
                continue

            for note in feature.qualifiers['note']:
                if note.startswith('tta leucine codon'):
                    tta_codons.append(feature)
                    break


        js_cluster = {}
        js_cluster['start'] = int(cluster.location.start) + 1
        js_cluster['end'] = int(cluster.location.end)
        js_cluster['idx'] = utils.get_cluster_number(cluster)
        js_cluster['orfs'] = convert_cds_features(record, features, annotations, options)
        js_cluster['borders'] = convert_cluster_border_features(borders)
        js_cluster['tta_codons'] = convert_tta_codons(tta_codons)
        js_cluster['type'] = utils.get_cluster_type(cluster)
        if 'probability' in cluster.qualifiers:
            js_cluster['probability'] = cluster.qualifiers['probability'][0]
        if options.input_type == 'prot':
            js_cluster['unordered'] = True
        js_cluster['knowncluster'] = "-"
        js_cluster['BGCid'] = "-"

        if 'knownclusterblast' in cluster.qualifiers:
            knownclusters = cluster.qualifiers['knownclusterblast']
            bestcluster = [kcluster for kcluster in knownclusters if kcluster.startswith('1.')]
            if not len(bestcluster) == 1:
                logging.warning("Error parsing best knowncluster hit; knownclusters array = %s. Possibly no significant hits to known biosynthetic gene clusters." % str(knownclusters))
            else:
                reObj = re.match('\d+\. (\S+)\t(.*)', bestcluster[0])
                js_cluster['knowncluster'] = reObj.group(2)
                js_cluster['BGCid'] = reObj.group(1)
                logging.debug('Found closest cluster "%s" for cluster no. %s' % (js_cluster['knowncluster'], utils.get_cluster_number(cluster)))
        js_clusters.append(js_cluster)

    return js_clusters

def convert_cds_features(record, features, annotations, options):
    """Convert CDS SeqFeatures to JSON"""
    js_orfs = []
    for feature in features:
        js_orf = {}
        js_orf['start'] = int(feature.location.start) + 1
        js_orf['end'] = int(feature.location.end)
        # Fix for files that have their coordinates the wrong way around
        if js_orf['start'] > js_orf['end']:
            js_orf['end'], js_orf['start'] = js_orf['start'], js_orf['end']
        js_orf['strand'] = feature.strand if feature.strand is not None else 1
        js_orf['locus_tag'] = utils.get_gene_id(feature)
        js_orf['type'] = get_biosynthetic_type(feature, annotations)
        js_orf['description'] = utils.ascii_string(get_description(record, feature, js_orf['type'], options))
        js_orfs.append(js_orf)
    return js_orfs


def convert_cluster_border_features(borders):
    js_borders = []
    for border in borders:
        border_note = border.qualifiers['note'][0]
        if not border_note.startswith('best prediction'):
            continue
        js_border = {}
        js_border['start'] = int(border.location.start) + 1
        js_border['end'] = int(border.location.end)
        # Always order the coordinates correctly
        if js_border['start'] > js_border['end']:
            js_border['end'], js_border['start'] = js_border['start'], js_border['end']
        js_border['tool'] = border.qualifiers['tool']
        js_borders.append(js_border)
    return js_borders


def convert_tta_codons(tta_codons):
    """Convert found TTA codon features to JSON"""
    js_codons = []
    for codon in tta_codons:
        js_codon = {}
        js_codon['start'] = int(codon.location.start) + 1
        js_codon['end'] = int(codon.location.end)
        js_codon['strand'] = codon.strand if codon.strand is not None else 1
        js_codons.append(js_codon)

    return js_codons

def get_description(record, feature, type_, options):
    "Get the description text of a feature"

    replacements = {
        'locus_tag': ", ".join(feature.qualifiers.get('locus_tag',['-'])),
        'protein_id': ", ".join(feature.qualifiers.get('protein_id',['-'])),
        'smcog': '-',
        'ecnumber': '-',
        'transport_blast_line': '',
        'smcog_tree_line': '',
        'searchgtr_line': '',
        'start': int(feature.location.start) + 1,
        'end': int(feature.location.end),
        'model_details': get_model_details(feature),
        'asf': ''
    }

    blastp_url = "http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&" \
                 "PROGRAM=blastp&BLAST_PROGRAMS=blastp&QUERY=%s&" \
                 "LINK_LOC=protein&PAGE_TYPE=BlastSearch"
    genomic_context_url = "http://www.ncbi.nlm.nih.gov/projects/sviewer/?" \
                          "Db=gene&DbFrom=protein&Cmd=Link&noslider=1&"\
                          "id=%s&from=%s&to=%s"
    template = '<span class="svgene-tooltip-bold">%(product)s</span><br>\n'
    template += 'Locus-tag: %(locus_tag)s; Protein-ID: %(protein_id)s<br>\n'
    if 'EC_number' in feature.qualifiers:
        template += "EC-number(s): %(ecnumber)s<br>\n"
    if options.smcogs:
        template += "smCOG: %(smcog)s<br>\n"
    if options.input_type == 'nucl':
        template += "Location: %(start)s - %(end)s<br><br>\n"
    if 'sec_met' in feature.qualifiers:
        template += '<span class="bold">Signature pHMM hits:</span><br>\n%(model_details)s<br>\n'

    if options.knownclusterblast:

        mibig_homology_path = glob(os.path.join(options.full_outputfolder_path, "knownclusterblast",
                                             "cluster*", utils.get_gene_acc(feature) + '_mibig_hits.txt'))
        if mibig_homology_path:
            mibig_homology_file = mibig_homology_path[0]
            generate_html_table(mibig_homology_file)
            html_file = mibig_homology_file.split('.txt')[0] + '.html'
            replacements['mibig_homology_path'] = html_file[len(options.full_outputfolder_path) + 1:]
            template += '<a href="%(mibig_homology_path)s" target="_new">MiBIG Hits</a><br><br>\n'
    template += """
%(transport_blast_line)s
%(searchgtr_line)s
<a href="%(blastp_url)s" target="_new">NCBI BlastP on this gene</a><br>
<a href="%(genomic_context_url)s" target="_new">View genomic context</a><br>
%(smcog_tree_line)s<br>"""
    if not get_ASF_predictions(feature) == "":
        template += '<span class="bold">Active Site Finder results:</span><br>\n%(asf)s<br><br>\n'
    template += """AA sequence: <a href="javascript:copyToClipboard('%(sequence)s')">Copy to clipboard</a><br>"""


    if not options.smcogs:
        del replacements['smcog']
    if options.input_type == 'prot':
        del replacements['start']
        del replacements['end']

    replacements['product'] = feature.qualifiers.get('product', ['-'])[0]
    if 'translation' in feature.qualifiers:
        sequence = feature.qualifiers['translation'][0]
    else:
        sequence = str(utils.get_aa_sequence(feature))
    replacements['blastp_url'] = blastp_url % sequence
    replacements['sequence'] = sequence
    if len(sequence) > 2000:
        len_seq = 30
    else:
        len_seq = (len(sequence) / 80) + 1
    replacements['len_seq'] = len_seq
    replacements['genomic_context_url'] = genomic_context_url % \
                    ( record.id,
                      max(feature.location.start - 9999, 0),
                      min(feature.location.end + 10000, len(record)) )
    if 'EC_number' in feature.qualifiers:
        replacements['ecnumber'] = ", ".join(feature.qualifiers.get('EC_number', ['-']))
    else:
        del replacements['ecnumber']

    if options.smcogs:
        for note in feature.qualifiers.get('note',[]):
            if note.startswith('smCOG:') and '(' in note:
                text = note[6:].split('(',1 )[0]
                smcog, desc = text.split(':', 1)
                desc = desc.replace('_', ' ')
                replacements['smcog'] = '%s (%s)' % (smcog, desc)
            elif note.startswith('smCOG tree PNG image:'):
                entry = '<a href="%s" target="_new">View smCOG seed phylogenetic tree with this gene</a>'
                url = note.split(':')[-1]
                replacements['smcog_tree_line'] = entry % url

    if type_ == 'transport':
        url = "http://blast.jcvi.org/er-blast/index.cgi?project=transporter;" \
              "program=blastp;database=pub/transporter.pep;" \
              "sequence=sequence%%0A%s" % sequence
        transport_blast_line = '<a href="%s" target="_new">TransportDB BLAST on this gene<br>' % url
        replacements['transport_blast_line'] = transport_blast_line

    if options.searchgtr_links.has_key(record.id + "_" + utils.get_gene_id(feature)):
        url = options.searchgtr_links[record.id + "_" + utils.get_gene_id(feature)]
        searchgtr_line = '<a href="%s" target="_new">SEARCHGTr on this gene<br>' % url
        replacements['searchgtr_line'] = searchgtr_line
    replacements['asf'] = get_ASF_predictions(feature)
    if replacements['asf'] == "":
        del replacements['asf']

    return template % replacements


def get_biosynthetic_type(feature, annotations):
    "Get the biosythetic type of a CDS feature"
    ann = 'other'
    for note in feature.qualifiers.get('note', []):
        if not note.startswith('smCOG:'):
            continue

        smcog = note[7:].split(':')[0]
        ann = annotations.get(smcog, 'other')

    if not 'sec_met' in feature.qualifiers:
        return ann

    for qual in feature.qualifiers['sec_met']:
        if not qual.startswith('Kind:'):
            continue
        return qual[6:]
    return ann

def get_model_details(feature):
    result = ""
    for note in feature.qualifiers.get('sec_met', []):
        if not note.startswith('Domains detected'):
            continue
        note = note[18:]
        result += note.replace(';', '<br>')

    return result

def get_ASF_predictions(feature):
    "check whether predictions from the active site finder module are annotated"


    ASFsec_met_quals = [sec_met_qual[16:] for sec_met_qual in feature.qualifiers.get('sec_met', [""]) if sec_met_qual.startswith("ASF-prediction")]
    result = "<br>\n".join(ASFsec_met_quals)

    return result

def load_cog_annotations():
    "Load the smCOG type annotations from a file"
    type_keys = {
        'B': 'biosynthetic-additional',
        'T': 'transport',
        'R': 'regulatory',
        'O': 'other'
    }
    annotations = {}
    for line in open(utils.get_full_path(__file__, 'cog_annotations.txt'), 'r'):
        line = line.strip()
        cog, _, type_ = line.split('\t', 3)
        annotations[cog] = type_keys.get(type_, 'other')

    return annotations
