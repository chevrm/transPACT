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
from pyquery import PyQuery as pq
import json
from os import path
import os
from antismash import utils
from antismash.config import get_config
from antismash.output_modules.html import js

def set_title(d, seq_id, num_clusters):
    d('title').text('%s - %s clusters - antiSMASH results' % (seq_id, num_clusters))

def set_version(d):
    d('#antismash-version').text(utils.get_version())


def set_colourscheme(d, options):
    """load the appropriate CSS file and logos"""
    d('link[rel=stylesheet]').attr.href = 'css/{}.css'.format(options.taxon)
    d('.antismash-logo').attr.src = 'images/{}_logo.png'.format(options.taxon)
    d('img[alt=home]').attr.src = 'images/{}_home.png'.format(options.taxon)
    d('img[alt=help]').attr.src = 'images/{}_help.png'.format(options.taxon)
    d('img[alt=about]').attr.src = 'images/{}_about.png'.format(options.taxon)
    d('img[alt=download]').attr.src = 'images/{}_download.png'.format(options.taxon)


def set_urls(d, options):
    if options.taxon == "fungi":
        base_url = options.urls.fungi_baseurl
    else:
        base_url = options.urls.bacteria_baseurl

    d('.main-link').attr.href = base_url
    d('.help-link').attr.href = base_url + '#!/help'
    d('.about-link').attr.href = base_url + '#!/about'


def set_download_links(d, seq_id, options):
    config = get_config()
    ul = d('#downloadoptions')

    item = pq('<li>')
    all_results = pq('<a>')
    all_results.text("Download all results")
    all_results.attr('href', '%s.zip' % seq_id)
    item.append(all_results)
    ul.append(item)

    if config.input_type == 'nucl':
        item = pq('<li>')
        xls = pq('<a>')
        xls.text("Download XLS overview file")
        xls.attr('href', '%s.geneclusters.xls' % seq_id)
        item.append(xls)
        ul.append(item)

        item = pq('<li>')
        embl = pq('<a>')
        embl.text("Download EMBL summary file")
        embl.attr('href', '%s.final.embl' % seq_id)
        item.append(embl)
        ul.append(item)

        item = pq('<li>')
        gbk = pq('<a>')
        gbk.text("Download GenBank summary file")
        gbk.attr('href', '%s.final.gbk' % seq_id)
        item.append(gbk)
        ul.append(item)

        if 'BiosynML' in options:
            item = pq('<li>')
            xml = pq('<a>')
            xml.text("Download BiosynML file")
            xml.attr('href', 'biosynML.xml')
            xml.attr('target', '_blank')
            item.append(xml)
            ul.append(item)

    else:
        item = pq('<li>')
        gbk = pq('<a>')
        gbk.text("Download GenPept summary file")
        gbk.attr('href', '%s.final.gp' % seq_id)
        item.append(gbk)
        ul.append(item)

    if 'logfile' in options:
        logfile_path = path.abspath(options.logfile)
        if path.dirname(logfile_path).startswith(options.full_outputfolder_path):
            rel_logfile_path = logfile_path[len(options.full_outputfolder_path)+1:]
            item = pq('<li>')
            logfile = pq('<a>')
            logfile.text("Download log file")
            logfile.attr('href', rel_logfile_path)
            item.append(logfile)
            ul.append(item)


def add_separator(d, seq_id, orig_id, options):
    overview_table = d('#cluster-overview>tbody')
    tr = pq('<tr>')
    tr.addClass('separator-row')
    # long description
    td = pq('<td>')
    td.addClass('separator-text')
    td.attr('colspan', '2')
    if options.input_type == 'nucl':
        id_text = seq_id
        if not orig_id == "":
            if len(orig_id)<40:
                id_text += " (original name was: %s)" % orig_id
            else:
                id_text += " (original name was: %s...)" % orig_id[:60]
        td.text('The following clusters are from record %s:' % id_text)
    else:
        td.text('The following clusters were found in your protein input:')
    tr.append(td)

    overview_table.append(tr)

def add_no_result_note(d, options):
    overview_table = d('#cluster-overview>tbody')
    tr = pq('<tr>')
    # clickable button
    td = pq('<td>')
    td.addClass('clbutton')
    td.addClass('separator')
    td.text('')
    tr.append(td)
    # long description
    td = pq('<td>')
    if options.input_type == 'nucl':
        td.text('No secondary metabolite clusters were found in the input sequence(s)')
    else:
        td.text('No secondary metabolite biosynthesis proteins were found in the input sequence(s)')
    tr.append(td)

    overview_table.append(tr)

def add_truncation_notice(d, options):
    header = d('#truncated')
    header.text(' (truncated to the first %s records)' % options.limit)

def add_cluster(d, cluster, seq_record, options, extra_data, odd, seq_id):
    add_cluster_button(d, cluster)
    add_overview_entry(d, cluster, odd)
    add_cluster_page(d, cluster, seq_record, options, extra_data, seq_id)

def add_cluster_button(d, cluster):
    li = pq('<li>')
    li.addClass('clbutton')
    li.addClass(cluster['type'])
    if cluster['type'].find('-') > -1:
        li.addClass('hybrid')
    li.addClass('cluster-%s' % cluster['idx'])
    a = pq('<a>')
    a.attr('href', '#cluster-%s' % cluster['idx'])
    a.text('%s' % cluster['idx'])
    li.append(a)
    # insert before the last element (the 'next' button)
    li.insert_before(d('#last-clbutton'))

def add_overview_entry(d, cluster, odd):
    overview_table = d('#cluster-overview>tbody')
    tr = pq('<tr>')
    if not odd:
        tr.addClass('even')
    # clickable button
    td = pq('<td>')
    td.addClass('clbutton')
    td.addClass(cluster['type'])
    if cluster['type'].find('-') > -1:
        td.addClass('hybrid')
    a = pq('<a>')
    a.attr('href', '#cluster-%s' % cluster['idx'])
    a.text('Cluster %s' % cluster['idx'])
    td.append(a)
    tr.append(td)
    # long description
    td = pq('<td>')
    subtypes = cluster['type'].split('-')
    i = 0
    for subtype in subtypes:
        a = pq('<a>')
        a.attr('href', "http://antismash.secondarymetabolites.org/help#{}".format(subtype))
        a.attr('target', '_blank')
        a.text(subtype.capitalize())
        td.append(a)
        if i < (len(subtypes) - 1):
            td.append('-')
        i += 1

    tr.append(td)
    # start
    td = pq('<td>')
    td.addClass('digits')
    td.text('%s' % cluster['start'])
    tr.append(td)
    # end
    td = pq('<td>')
    td.addClass('digits')
    td.text('%s' % cluster['end'])
    tr.append(td)
    # closest cluster match BGCid description
    td = pq('<td>')
    td.text(cluster['knowncluster'])
    tr.append(td)

    td = pq('<td>')
    if cluster['BGCid'] != '-':
        a = pq('<a>')
        a.attr('href', "http://mibig.secondarymetabolites.org/repository/{}/index.html".format(cluster['BGCid'].split('_')[0]))
        a.attr('target', '_blank')
        a.text(cluster['BGCid'])
        td.append(a)
    else:
       td.text(cluster['BGCid'])
    tr.append(td)


    overview_table.append(tr)

def add_cluster_page(d, cluster, seq_record, options, extra_data, seq_id):

    handlers = find_plugins_for_cluster(options.plugins, cluster)

    cluster_rec = utils.get_cluster_by_nr(seq_record, cluster['idx'])

    rules = get_detection_rules(cluster_rec)

    page = pq('<div>')
    page.addClass('page')
    page.attr('id', 'cluster-%s' % cluster['idx'])
    header = pq('<h3>')
    header.text('%s - Cluster %s - %s' % (seq_record.name, cluster['idx'], cluster['type'].capitalize()))
    page.append(header)

    sidepanel = None
    for handler in handlers:
        sidepanel = handler.generate_sidepanel(cluster, seq_record, options, sidepanel)

    if sidepanel is not None:
        page.append(sidepanel)

    content = pq('<div>')
    content.addClass('content')

    description = pq('<div>')
    description.addClass('description-container')
    desc_header = pq('<h3>')
    desc_header.text('Gene cluster description')
    description.append(desc_header)

    cluster_download = pq('<div>')
    cluster_download.addClass('cluster-download')
    description.append(cluster_download)
    dl_link = pq('<a>')
    dl_link.attr('href', '%s.cluster%03d.gbk' % (seq_id, cluster['idx']))
    dl_link.text('Download cluster GenBank file')
    cluster_download.append(dl_link)

    desc_text = pq('<div>')
    desc_text.addClass('description-text')
    if options.input_type == 'nucl':
        text = seq_record.name +' - Gene Cluster %(idx)s. Type = %(type)s. Location: %(start)s - %(end)s nt. '
    else:
        text = seq_record.name + '- Gene Cluster %(idx)s. Type = %(type)s. '
    if 'probability' in cluster:
        text += 'ClusterFinder probability: %(probability)s. '
    text += 'Click on genes for more information.'
    desc_text.text(text % cluster)
    description.append(desc_text)
    rules_header = pq('<a>')
    rules_header.addClass('cluster-rules-header')
    rules_header.attr('id', 'cluster-%s-rules-header'% cluster['idx'])
    rules_header.attr('href', '#cluster-%s' % cluster['idx'])
    rules_header.text('Show pHMM detection rules used')
    description.append(rules_header)

    detection_rules = pq('<div>')
    detection_rules.addClass('cluster-rules')
    detection_rules.attr('id', 'cluster-%s-rules' % cluster['idx'])
    detection_rules.html('<br>'.join(rules))
    description.append(detection_rules)


    desc_svg = pq('<div>')
    desc_svg.attr('id', 'cluster-%s-svg' % cluster['idx'])
    description.append(desc_svg)

    content.append(description)

    if options.input_type == 'nucl':
        legend = pq('<div>')
        legend.addClass('legend')
        legend_header = pq('<h4>')
        legend_header.text('Legend:')
        legend.append(legend_header)
        legend_text = pq('<div>')
        if not options.smcogs:
            legend_text.append("Only available when smCOG analysis was run")
        legend_text.append(generate_legend_entry('legend-type-biosynthetic', 'core biosynthetic genes'))
        legend_text.append(generate_legend_entry('legend-type-biosynthetic-additional', 'additional biosynthetic genes'))
        legend_text.append(generate_legend_entry('legend-type-transport', 'transport-related genes'))
        legend_text.append(generate_legend_entry('legend-type-regulatory', 'regulatory genes'))
        legend_text.append(generate_legend_entry('legend-type-other', 'other genes'))
        if options.tta:
            legend_text.append(generate_legend_entry('legend-tta-codon', 'TTA codon'))
        if options.cassis:
            legend_text.append(generate_legend_entry('legend-border-cassis', 'cluster extent as predicted by CASSIS'))
        if options.borderpredict:
            legend_text.append(generate_legend_entry('legend-border-clusterfinder', 'cluster extent as predicted by ClusterFinder'))
        legend.append(legend_text)
        content.append(legend)

    details = None
    da = None
    for handler in handlers:
        if "generate_details_div" in dir(handler):
            details = handler.generate_details_div(cluster, seq_record, options,
                                                   extra_data['js_domains'], details)
        if 'generate_domain_alignment_div' in dir(handler) and options.transatpks_da:
            da = handler.generate_domain_alignment_div(cluster, seq_record, options, da)

    if details is not None:
        content.append(details)

    if da is not None:
        content.append(da)

    if options.clusterblast:
        top_ten_clusters = cluster_rec.qualifiers.get('clusterblast', [])

        cb = pq('<div>')
        cb.addClass('clusterblast')
        cb_header = pq('<h3>')
        cb_header.text("Homologous gene clusters")
        cb.append(cb_header)
        cb_control = pq('<div>')
        cb.append(cb_control)
        if len(top_ten_clusters) == 0:
            cb_download = pq('No significant ClusterBlast hits found.')
            cb_control.append(cb_download)
        else:
            cb_select = pq('<select>')
            cb_select.attr('id', 'clusterblast-%s-select' % cluster['idx'])
            cb_select.addClass('clusterblast-selector')
            cb_control.append(cb_select)
            opt = pq('<option>')
            opt.attr('value', path.join('svg', 'clusterblast%s_all.svg' % cluster['idx']))
            opt.text('All hits')
            cb_select.append(opt)
            for i in range(1, options.nclusters+1):
                svg_file = path.join('svg', 'clusterblast%s_%s.svg' % (cluster['idx'], i))
                full_path = path.join(options.outputfoldername, svg_file)
                if path.exists(full_path):
                    opt = pq('<option>')
                    opt.attr('value', svg_file)
                    opt_text = 'Cluster %s hit %s' % (cluster['idx'], i)
                    if len(top_ten_clusters) >= i:
                        opt_text = top_ten_clusters[i-1].split('\t')[1]
                    opt.text(opt_text)
                    cb_select.append(opt)
                else:
                    logging.debug("failed to find %r" % full_path)
            cb_download = pq('<button>');
            cb_download.attr('id', 'clusterblast-%s-download' % cluster['idx'])
            cb_download.text('Download graphic')
            cb_control.append(cb_download)

        cb_svg = pq('<div>')
        cb_svg.attr('id', 'clusterblast-%s-svg'% cluster['idx'])
        cb.append(cb_svg)
        content.append(cb)

    if options.subclusterblast:
        top_ten_clusters = cluster_rec.qualifiers.get('subclusterblast', [])

        cb = pq('<div>')
        cb.addClass('subclusterblast')
        cb_header = pq('<h3>')
        cb_header.text("Homologous subclusters")
        cb.append(cb_header)
        cb_control = pq('<div>')
        cb.append(cb_control)
        cb_select = pq('<select>')
        cb_select.attr('id', 'subclusterblast-%s-select' % cluster['idx'])
        cb_select.addClass('clusterblast-selector')
        cb_control.append(cb_select)
        opt = pq('<option>')
        opt.attr('value', path.join('svg', 'subclusterblast%s_all.svg' % cluster['idx']))
        opt.text('All hits')
        cb_select.append(opt)
        subclusters_added = 0
        for i in range(1, options.nclusters+1):
            svg_file = path.join('svg', 'subclusterblast%s_%s.svg' % (cluster['idx'], i))
            full_path = path.join(options.outputfoldername, svg_file)
            if path.exists(full_path):
                opt = pq('<option>')
                opt.attr('value', svg_file)
                opt_text = 'Cluster %s hit %s' % (cluster['idx'], i)
                if len(top_ten_clusters) >= i:
                    opt_text = top_ten_clusters[i-1].split('\t')[1].replace('_', ' ')
                opt.text(opt_text)
                cb_select.append(opt)
                subclusters_added += 1
            else:
                logging.debug("failed to find %r" % full_path)

        cb_svg = pq('<div>')
        cb_svg.attr('id', 'subclusterblast-%s-svg'% cluster['idx'])
        cb.append(cb_svg)
        if path.exists(path.join(options.outputfoldername, 'svg',
                       'subclusterblast%s_all.svg' % cluster['idx'])) and \
           subclusters_added > 0:
            cb_download = pq('<button>');
            cb_download.attr('id', 'subclusterblast-%s-download' % cluster['idx'])
            cb_download.text('Download graphic')
            cb_control.append(cb_download)
            content.append(cb)

    if options.knownclusterblast:
        top_ten_clusters = cluster_rec.qualifiers.get('knownclusterblast', [])

        cb = pq('<div>')
        cb.addClass('knownclusterblast')
        cb_header = pq('<h3>')
        cb_header.text("Homologous known gene clusters")
        cb.append(cb_header)
        cb_control = pq('<div>')
        cb.append(cb_control)
        cb_select = pq('<select>')
        cb_select.attr('id', 'knownclusterblast-%s-select' % cluster['idx'])
        cb_select.addClass('clusterblast-selector')
        cb_control.append(cb_select)
        opt = pq('<option>')
        opt.attr('value', path.join('svg', 'knownclusterblast%s_all.svg' % cluster['idx']))
        opt.text('All hits')
        cb_select.append(opt)
        knownclusters_added = 0
        for i in range(1, options.nclusters+1):
            svg_file = path.join('svg', 'knownclusterblast%s_%s.svg' % (cluster['idx'], i))
            full_path = path.join(options.outputfoldername, svg_file)
            if path.exists(full_path):
                opt = pq('<option>')
                opt.attr('value', svg_file)
                opt_text = 'Cluster %s hit %s' % (cluster['idx'], i)
                if len(top_ten_clusters) >= i:
                    opt_text = top_ten_clusters[i-1].split('\t')[1].replace('_', ' ')
                opt.text(opt_text)
                cb_select.append(opt)
                knownclusters_added += 1
            else:
                logging.debug("failed to find %r" % full_path)

        cb_svg = pq('<div>')
        cb_svg.attr('id', 'knownclusterblast-%s-svg'% cluster['idx'])
        cb.append(cb_svg)
        if path.exists(path.join(options.outputfoldername, 'svg',
                       'knownclusterblast%s_all.svg' % cluster['idx'])) and \
           knownclusters_added > 0:
            cb_download = pq('<button>');
            cb_download.attr('id', 'knownclusterblast-%s-download' % cluster['idx'])
            cb_download.text('Download graphic')
            cb_control.append(cb_download)
            content.append(cb)

    page.append(content)
    d('.page:last').after(page)


def write_geneclusters_js(records, output_dir, extra_data):
    with open(path.join(output_dir, 'geneclusters.js'), 'w') as h:
        geneclusters = {}
        for record in records:
            for cluster in record['clusters']:
                idx = cluster['idx']
                geneclusters['cluster-%s' % idx] = cluster
        h.write('var geneclusters = %s;\n' %
                json.dumps(geneclusters, indent=4))

        js_domains = {}
        for domain in extra_data['js_domains']:
            idx = domain['id'].split('-')[1]
            js_domains['cluster-%s' % idx] = domain
        h.write('var details_data = %s;\n' %
                json.dumps(js_domains, indent=4))

def generate_webpage(seq_records, options):
    d = pq(filename=utils.get_full_path(__file__, 'index.tpl'), parser='html')

    num = count_all_clusters(seq_records)
    set_title(d, seq_records[0].id, num)
    set_colourscheme(d, options)
    set_urls(d, options)
    set_version(d)
    set_download_links(d, seq_records[0].id, options)

    generate_searchgtr_htmls(seq_records, options)
    records = js.convert_records(seq_records, options)

    extra_data = dict(js_domains=[], clusterblast_clusters=[],
                      subclusterblast_clusters=[], knownclusterblast_clusters=[])

    if 'triggered_limit' in options and options.triggered_limit:
        add_truncation_notice(d, options)

    records_written = 0

    for i in range(len(records)):
        odd = True
        records[i]['seq_id'] = utils.ascii_string(records[i]['seq_id'])
        if len(records[i]['clusters']) > 0:
            add_separator(d, records[i]['seq_id'], records[i]['orig_id'], options)
        for cluster in records[i]['clusters']:
            add_cluster(d, cluster, seq_records[i], options, extra_data, odd, seq_records[0].id)
            records_written += 1
            odd = not odd

    if records_written == 0:
        add_no_result_note(d, options)

    write_geneclusters_js(records, options.outputfoldername, extra_data)

    with open(path.join(options.outputfoldername, 'index.html'), 'w') as h:
        h.write('<!doctype html>\n')
        h.write(d.outerHtml())

def find_plugins_for_cluster(plugins, cluster):
    "Find a specific plugin responsible for a given gene cluster type"
    product = cluster['type']
    handlers = []
    for plugin in plugins:
        if plugin.will_handle(product):
            handlers.append(plugin)

    return handlers

def count_all_clusters(records):
    "Count the number of cluster across all records"
    num = 0
    for record in records:
        num += len(utils.get_cluster_features(record))

    return num

def generate_legend_entry(css_class, description):
    "Create an entry for a legend, based on a css clas of the svg element"
    entry = pq('<div>')
    box = pq('<div>')
    box.addClass('legend-field')
    box.addClass(css_class)
    entry.append(box)

    label = pq('<div>')
    label.addClass('legend-label')
    label.text(description)
    entry.append(label)

    return entry

def get_detection_rules(cluster_rec):
    rules = []
    for note in cluster_rec.qualifiers['note']:
        if note.startswith("Detection rule(s)"):
            rules.extend([rule.strip().replace('&', '&amp;') for rule in note[41:].split(';')])

    return rules

def load_searchgtr_search_form_template():
    #Create folder for SEARCHGTR HTML files, load search form template
    searchgtrformtemplate = open(path.join(utils.get_full_path(__file__, ''), "searchgtr_form.html"),"r")
    searchgtrformtemplate = searchgtrformtemplate.read()
    searchgtrformtemplate = searchgtrformtemplate.replace("\r","\n")
    searchgtrformtemplateparts = searchgtrformtemplate.split("FASTASEQUENCE")
    return searchgtrformtemplateparts

def generate_searchgtr_htmls(seq_records, options):
    #Generate lists of COGs that are glycosyltransferases or transporters
    gtrcoglist = ['SMCOG1045','SMCOG1062','SMCOG1102']
    searchgtrformtemplateparts = load_searchgtr_search_form_template()
    options.searchgtr_links = {}
    for seq_record in seq_records:
        smcogdict, _ = utils.get_smcog_annotations(seq_record)
        for feature in utils.get_cds_features(seq_record):
            gene_id = utils.get_gene_id(feature)
            if smcogdict.has_key(gene_id):
                smcog = smcogdict[gene_id]
                if smcog in gtrcoglist:
                
                    if not os.path.exists(options.full_outputfolder_path + os.sep + "html"):
                        os.mkdir(options.full_outputfolder_path + os.sep + "html")
                    formfileloc = options.full_outputfolder_path + os.sep + "html" + os.sep + utils.get_gene_id(feature) + "_searchgtr.html"
                    link_loc = "html" + os.sep + utils.get_gene_id(feature) + "_searchgtr.html"
                    options.searchgtr_links[seq_record.id + "_" + gene_id] = link_loc
                    formfile = open(formfileloc, "w")
                    specificformtemplate = searchgtrformtemplateparts[0].replace("GlycTr",gene_id)
                    formfile.write(specificformtemplate)
                    formfile.write("%s\n%s" % (gene_id, utils.get_aa_sequence(feature)))
                    formfile.write(searchgtrformtemplateparts[1])
                    formfile.close()
