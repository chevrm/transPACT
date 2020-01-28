import re
import collections
import logging
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import xml.etree.ElementTree as ET
ET.register_namespace("", "http://www.w3.org/2000/svg")
ET.register_namespace('xlink', "http://www.w3.org/1999/xlink")
import Longest_Common_Subsequence as LCS
import get_speciesName_in_figures as GSN
import get_module_contig_info as MCI

from antismash.utils import sort_XonY

def modify_svg(f, node_tooltip_description, link_tooltip_description, f_out, cluster_id):
    t = ET.parse(f)
    root = t.getroot()
    #change size of figure
    root.attrib['width'] = '1440pt'
    root.attrib['height'] = '600pt'
    root.attrib['viewBox'] = '500 300 2880 600'
    #get all elements which are either for links and nodes"
    elem_handel_link = []
    for elem in root.iter('{http://www.w3.org/2000/svg}g'):
        if 'id' in elem.attrib.keys():
            if "LineCollection_" in elem.get('id'):
                elem_handel_link.append(elem)
    elem_handel_node = []
    pos_node = []
    for elem in root.iter('{http://www.w3.org/2000/svg}g'):
        if 'id' in elem.attrib.keys():
            # PathCollection_1 is the * to denote the query
            if "PathCollection_" in elem.get('id') and elem.get('id') != "PathCollection_1":
                node_use = [u for u in elem.iter() if u.tag == "{http://www.w3.org/2000/svg}use"]
                if not node_use:
                    continue
                position1 = float(node_use[0].get('y'))
                if len(node_use) < 2:
                    elem_handel_node.append(elem)
                    pos_node.append(position1)
                    continue
                position2 = float(node_use[1].get('y'))
                # this is only true for the legend
                if position1 == position2:
                    elem_handel_node.append(elem)
                    pos_node.append(position1)
    elem_handel_node_order = sort_XonY(x=elem_handel_node, y=pos_node, reverse=True)
    id_start = 0
    link_start = 0
    for elem in elem_handel_link:
        if elem.get('id') == "LineCollection_1":
            lines = [l for l in elem.iter() if l.tag == "{http://www.w3.org/2000/svg}path"]
            for p in lines:
                p.attrib['class'] = 'domainalign-orf'
                p.attrib['locus_tag'] = ''
                p.attrib['description'] = link_tooltip_description[link_start]
                p.attrib['id'] = '-'.join(['domainalign', cluster_id, str(id_start)])
                link_start += 1
                id_start += 1
    node_start = 0
    for elem in elem_handel_node_order:
        if 'PathCollection_' in elem.get('id'):
            nodes = [n for n in elem.iter() if n.tag == "{http://www.w3.org/2000/svg}g" and n != elem]
            for n in nodes:
                n.attrib['class'] = 'domainalign-orf'
                n.attrib['locus_tag'] = ''
                n.attrib['description'] = node_tooltip_description[node_start]
                n.attrib['id'] = '-'.join(['domainalign', cluster_id, str(id_start)])
                node_start += 1
                id_start += 1
    t.write(f_out)

#-------------------------------------------------------------------------------------------------
#-- function: plot_dendrogram_upgma
#-- This is the function to genetate a hierarchical clustering plot
#--
#-- The input file of this function is a distance matrix, where each entry is the similiary score
#-- between two assembly lines. The similarity score is determined by considering:
#-- i) the composition of substarte-specificity of the domains
#-- ii) the  order of KS domains
#-- iii) sequence similarity between two KS domains
#--------------------------------------------------------------------------------------------------
def plot_dendrogram_upgma(dist_df, outfile):
    labelList = list(dist_df.dtypes.index)
    Linkage_matrix = sch.linkage(dist_df, method = 'average')
    fig_h = dist_df.shape[0]
    fig = plt.figure()
    fig.set_size_inches(w = 3, h = 0.8*fig_h)
    plt.title('UPGMA clustering')
    dn = sch.dendrogram(
        Linkage_matrix,
        orientation = 'left',
        labels = labelList,
        get_leaves = True,
        no_labels=True
    )
    leaves = dn['ivl']
    leaves_reformat =  [i.split('(')[0] for i in leaves]
    plt.axis('off')
    plt.savefig(outfile, format = 'svg')
    return leaves_reformat

def get_bgc_dendrogram_upgma(dist_df):
    labelList = list(dist_df.dtypes.index)
    Linkage_matrix = sch.linkage(dist_df, method = 'average')
    dn = sch.dendrogram(
        Linkage_matrix,
        orientation = 'left',
        labels = labelList,
        get_leaves = True,
        no_labels=True,
        no_plot = True
    )
    leaves = dn['ivl']
    leaves_reformat =  [i.split('(')[0] for i in leaves]
    return leaves_reformat

def extract_KS_anno_unique(id_list, new_bgc_id, new_cluster, new_cluster_index, annoKS, annoKS_n, KSnames):
    id_old_list = [id for id in id_list if id not in new_bgc_id]
    KS_id = annoKS_n[0]
    output_anno = {}
    output_index = {}

    for i in id_old_list:
        compd_id = i
        compd_i = KSnames[compd_id]

        compd_ks_index = zip(compd_i.split('\t')[1:], KS_id.split('\t')[1:])
        compd_ks = [x for x,y in compd_ks_index if x != 'NA' ][1:]
        compd_index = [y for x,y in compd_ks_index if x != 'NA'][1:]
        output_anno[compd_id] = ';'.join(compd_ks)
        output_index[compd_id] = compd_index
    #-- the assembly lines are filtered only based on the domian composition, the domain sequence similarity are not considered
    output_anno_val = list(set(output_anno.values()))
    output_anno_reverse = {}
    for anno in output_anno_val:
        ID = []
        for x,y in output_anno.items():
            if y == anno:
                ID.append(x)
        output_anno_reverse[anno] = ID
    output_anno_mergeID = {}
    repeatID_dict = {}
    output_index_mergeID = {}
    for k in output_anno_reverse.keys():
        val = output_anno_reverse[k]
        if len(val) == 1:
            new_key = val[0]
            output_anno_mergeID[new_key] = k.split(';')
            output_index_mergeID[new_key] = output_index[new_key]
        else:
            new_key = ''.join([val[0], '(extra ', str(len(val) -1), ')'])
            repeatID_dict[new_key] = ';'.join(val)
            output_anno_mergeID[new_key] = k.split(';')
            output_index_mergeID[new_key] = output_index[val[0]]
    if len(id_old_list) < len(id_list):
        new_cluster_dict = {}
        new_cluster_index_dict = {}
        for cluster in new_cluster:
            k = cluster.split('\n')[0]
            v = cluster.split('\n')[1].split(',')
            new_cluster_dict[k] = v
        for cluster in new_cluster_index:
            k = cluster.split('\n')[0]
            v = cluster.split('\n')[1].split(',')
            new_cluster_index_dict[k] = v
        for k_id in new_bgc_id:
            output_anno_mergeID[k_id] = new_cluster_dict[k_id]
            output_index_mergeID[k_id] = new_cluster_index_dict[k_id]
    return output_anno_mergeID, output_index_mergeID

def create_node(input, KSindex, KSanno, dist):
    '''
    :param input:
    :param KSindex:
    :param KSanno:
    :param dist:
    :return:
    KSanno = KS_anno
    KSindex = KS_index
    input = novel[0:5]
    dist = 0.3
    nodes, labels, positions, xlim, ylim, clade = create_node(input = novel[0:2], KSindex = KS_index, KSanno = KS_anno, dist = 0.3)
    '''
    node_anno = []
    node_labels = []
    node_y_loc = []
    node_x_loc = []
    node_len = []
    node_clade = []
    legend_x_loc = {}
    legend_y_loc = {}
    node_index_dict = {} # this is the dictionary, where key is the compound id, and node is the node index, {'dorrigocin_migrastatin|GQ274953': [0,1,2]}
    start = 1
    node_start = 0
    compd_node_label = {}
    for i in input:
        label_i = KSanno[i]
        label_i_adj = []
        for l in label_i:
            if 'noInfo' in l or 'not_conserved' in l or 'clade_not_conserved' in l:
                label_i_adj.append('NA')
            else:
                label_i_adj.append(l)
        anno_i = ["|".join((i, v)) for v in KSindex[i]]
        node_anno = node_anno + anno_i
        node_labels = node_labels + label_i_adj
        node_clade = node_clade + KSanno[i]
        l = len(KSindex[i])
        node_index = range(node_start, node_start+l)
        node_index_dict[i] = node_index
        node_len.append(l)
        node_y_loc = node_y_loc + [start]*l
        node_x_loc = node_x_loc + [1+dist*n for n in range(0,l)]
        compd_node_label[i] = node_index
        start += line_gap
        node_start = max(node_index)+1
        legend_x_loc[i] = 1
        legend_y_loc[i] = start-line_gap
    node_loc = zip(node_x_loc, node_y_loc)
    node = range(0, sum(node_len))
    node_loc_dict = dict(list(enumerate(node_loc)))
    node_labels_dict = dict(list(enumerate(node_labels)))
    node_clade_dict = dict(list(enumerate(node_clade)))
    node_anno_dict = dict(list(enumerate(node_anno)))
    clade_output = {}
    for i in input:
        index = compd_node_label[i]
        dict_i = {k:node_clade_dict[k] for k in index}
        clade_output[i] = dict_i
    nr_unique_phyloClade = len([l for l in list(set([v for sublist in KSanno.values() for v in sublist])) if 'noInfo' not in l])
    xlim_a = min(node_x_loc)-3
    xlim_b = max(node_x_loc)+0.5
    ylim_a = min(node_y_loc)-1
    ylim_b = max(node_y_loc)+1
    if round(float(nr_unique_phyloClade)/(2*ylim_b-1)) > float(nr_unique_phyloClade)/(2*ylim_b-1):
        xlim_b_mod = xlim_b + 2 * round(float(nr_unique_phyloClade)/(2*ylim_b-1))
    else:
        xlim_b_mod = xlim_b + 2 * (round(float(nr_unique_phyloClade)/(2*ylim_b-1)) + 1)+1
    return node, node_labels_dict, node_loc_dict, [xlim_a, xlim_b_mod], [ylim_a, ylim_b], clade_output, legend_x_loc, legend_y_loc, [xlim_b, ylim_b-1], node_index_dict, node_anno_dict

def color_nodes_phyloClade(KSanno, color_code_file):
    #-- create a dict of color codes, {'phyloClade_name': '#4169E1', ....}
    color_code = [c for c in open(color_code_file, 'r').read().split('\n')[1:] if c != '']
    color_code_dict = {}
    for c in color_code:
        k = c.split('\t')[1]
        v = c.split('\t')[4]
        color_code_dict[k] = v
    #-- create the dict containing color code for nodes in each assembly line, {'assemly_line_ID1': ['#4169E1', '#00008B', ...], ''}
    KS_color = {}
    node_edge_color = {}
    for k in KSanno.keys():
        compd_phyloClade = KSanno[k]
        compd_color = []
        compd_edge_color = []
        for p in compd_phyloClade:
            if p not in color_code_dict.keys():
                compd_color.append('#FFFFFF')
                compd_edge_color.append('#000000')
            else:
                color = color_code_dict[p]
                compd_color.append(color)
                compd_edge_color.append(color)
        KS_color[k] = compd_color
        node_edge_color[k] = compd_edge_color
    return KS_color, node_edge_color

def get_edge_pairwise(compd1, compd2, clade_node_label, KSanno):
    '''
    :param compd1:
    :param compd2:
    :param clade_node_label:
    :param KSanno:
    :return:
    compd1 = novel[0]
    compd2 = novel[1]
    clade_node_label = clade
    KSanno = KS_anno
    edge_list = get_edge_pairwise(compd1 = novel[0], compd2 = novel[1], clade_node_label = clade, KSanno = KS_anno)
    '''
    compd1_KS = KSanno[compd1]
    compd2_KS = KSanno[compd2]
    compd1_node = clade_node_label[compd1]
    compd2_node = clade_node_label[compd2]
    #we attached a number for 'noInfo', so that no match will be count for noInfo-noInfo pairs
    affli = 0
    for i in range(0, len(compd1_KS)):
        if compd1_KS[i] == 'noInfo':
            compd1_KS[i] = 'noInfo'+ str(affli)
            affli += 1
    for i in range(0, len(compd2_KS)):
        if compd2_KS[i] == 'noInfo':
            compd2_KS[i] = 'noInfo' + str(affli)
            affli += 1
    a, b = LCS.alignSequences(compd1_KS, compd2_KS)
    pointer_compd1 = []
    pointer_compd2 = []
    compd1_index = sorted(compd1_node.keys())
    compd2_index = sorted(compd2_node.keys())
    for i in range(0,len(a)):
        if a[i] != '' and b[i] != '':
            p1 = i - len([m for m in range(0,i) if a[m] == ''])
            p2 = i - len([n for n in range(0,i) if b[n] == ''])
            pointer_compd1.append(p1)
            pointer_compd2.append(p2)
    compd1_node_aligned = [compd1_index[index] for index in pointer_compd1]
    compd2_node_aligned = [compd2_index[index] for index in pointer_compd2]
    edge = zip(compd1_node_aligned, compd2_node_aligned)
    return edge

def get_node_edge_multiAlignment(seq_ordered, KS_anno, KS_index, dist):
    nodes, labels, positions, xlim, ylim, clade, legend_x, legend_y, pos_node, node_index_dict, node_anno_dict = create_node(input = seq_ordered, KSindex = KS_index, KSanno = KS_anno, dist = dist)
    seq_pairs = zip(range(0,len(seq_ordered)-1), range(1,len(seq_ordered)))
    edge_output = []
    for x,y in seq_pairs:
        edge_i = get_edge_pairwise(compd1 = seq_ordered[x], compd2 = seq_ordered[y], clade_node_label = clade, KSanno = KS_anno)
        edge_output = edge_output + edge_i
    return nodes, labels, positions, xlim, ylim, edge_output, legend_x, legend_y, pos_node, node_index_dict, node_anno_dict

def get_legend(KSanno, color_code_file, pos_node, nodes_max):
    #-- create a dict of color code and its corresponding function, {'phyloClade_name': ['#4169E1', 'glycine'], ....}
    color_code = [c for c in open(color_code_file, 'r').read().split('\n')[1:] if c != '']
    color_fun_dict = {}
    for c in color_code:
        k = c.split('\t')[1]
        v1 = c.split('\t')[4]
        v2 = c.split('\t')[3]
        color_fun_dict[k] = [v1, v2]
    unique_phyloClade = extract_unique_phylo_clade(KSanno)
    cleanup_unique_phylo_clade(unique_phyloClade)
    legend_annotate = []
    legend_color = []
    legend_nodes_index = []
    legend_x_pos = []
    legend_y_pos = []
    start = 0
    num = int((pos_node[1] - 1)/line_gap + 1)
    estimate_col = float(len(unique_phyloClade)) / (2 * num - 1)
    if round(estimate_col) >= estimate_col:
        col_nr = int(round(estimate_col))
    else:
        col_nr = int(round(estimate_col) + 1)
    for k in range(col_nr):
        x_pos = pos_node[0] + 2 * k
        for j in range(2 * num - 1):
            index = j + k * (2 * num - 1)
            if index <= len(unique_phyloClade)-1:
                clade_name = unique_phyloClade[index]
                y_pos = pos_node[1] - (line_gap*0.7) * j
                legend_annotate.append(color_fun_dict[clade_name][1])
                legend_color.append(color_fun_dict[clade_name][0])
                legend_x_pos.append(x_pos)
                legend_y_pos.append(y_pos)
                legend_nodes_index.append(start+nodes_max)
                start += 1
            else:
                break
    nodes_legend = legend_nodes_index
    nodes_legend_labels = dict(zip(legend_nodes_index, unique_phyloClade))
    nodes_legend_pos = dict(zip(legend_nodes_index, zip(legend_x_pos, legend_y_pos)))
    nodes_legend_color = dict(zip(legend_nodes_index, legend_color))
    nodes_legend_annotation = dict(zip(legend_nodes_index, legend_annotate))
    return nodes_legend, nodes_legend_labels, nodes_legend_pos, nodes_legend_color, nodes_legend_annotation

def extract_unique_phylo_clade(ks_anno):
    # get unique phylogeney clades in the extracted assembly lines
    clade_set = set()
    annotation_lists = ks_anno.values()
    for inner_list in annotation_lists:
        for val in inner_list:
            if 'noInfo' in val:
                continue
            if 'tree_not_conserved' in val:
                continue
            if 'clade_not_conserved' in val:
                continue
            if 'out_group_like' in val:
                continue
            clade_set.add(val)
    unique_phylo_clade = list(clade_set)
    unique_phylo_clade.sort()
    return unique_phylo_clade

def cleanup_unique_phylo_clade(unique_phylo_clade):
    try:
        # We shouldn't need a loop, but better safe than sorry
        idx = unique_phylo_clade.index('II_c5_g1')
        while idx > -1:
            unique_phylo_clade[idx] = 'II_c5'
            idx = unique_phylo_clade.index('II_c5_g1', idx + 1)
    except ValueError:
        pass
    try:
        # see above, better safe than sorry
        idx = unique_phylo_clade.index('II_c5_g2')
        while idx > -1:
            unique_phylo_clade[idx] = 'II_c5'
            idx = unique_phylo_clade.index('II_c5_g2', idx + 1)
    except ValueError:
        pass
    try:
        first_idx = unique_phylo_clade.index('II_c5')
        if first_idx > -1:
            # we want to keep the first one
            idx = unique_phylo_clade.index('II_c5', first_idx + 1)
            while idx > -1:
                del unique_phylo_clade[idx]
                idx = unique_phylo_clade.index('II_c5', first_idx + 1)
    except ValueError:
        pass

def get_figureInfo_domain_loc_edge(node_index_dict, positions, domain_loc_edge, KSindex):
    loc_output = []
    edge_left_output = []
    edge_right_output = []
    for k in domain_loc_edge.keys():
        KS_id = KSindex[k]
        node_id = node_index_dict[k]
        KS_node_id = dict(zip(KS_id, node_id))
        loc = domain_loc_edge[k][0]
        edge_left = domain_loc_edge[k][1]
        edge_right = domain_loc_edge[k][2]
        loc_T = []
        for x,y in zip(node_id, loc):
            if y == 'T':
                loc_T.append(x)
        loc_T_plus = [i+1 for i in loc_T]
        loc_output_k = zip(loc_T, loc_T_plus)
        loc_output = loc_output + loc_output_k
        if edge_left != []:
            for i in edge_left:
                edge_left_output.append(KS_node_id[i])
        if edge_right != []:
            for i in edge_right:
                edge_right_output.append(KS_node_id[i])
    edge_sign_axisX = []
    edge_sign_axisY = []
    if edge_right_output != []:
        for r in edge_right_output:
            X = positions[r][0] + 0.15
            Y = positions[r][1]
            edge_sign_axisX.append(X)
            edge_sign_axisY.append(Y)
    if edge_left_output != []:
        for l in edge_left_output:
            X = positions[l][0] - 0.15
            Y = positions[l][1]
            edge_sign_axisX.append(X)
            edge_sign_axisY.append(Y)
    return loc_output, edge_sign_axisX, edge_sign_axisY

def get_tooltip_links_nodes(node_dict, edge_list, seq_simScore):
    DM_input = [d for d in open(seq_simScore, 'r').read().split('\n') if d != '']
    index = [0, 1, 9]
    domain_name = []
    score_len = []
    for l in DM_input:
        s = len(l.split(','))
        score_len.append(s)
    for l in DM_input:
        domain = "|".join([l.split(",")[0].split("|")[i] for i in index])
        domain_name.append(domain)
    domain_index = collections.defaultdict(list)
    for n, i in zip(domain_name, range(0, len(domain_name))):
        domain_index[n].append(i)
    node_tooltip_description = []
    for i in range(len(node_dict)):
        node_id = node_dict[i]
        row_ix = domain_index[node_id][0]
        info = DM_input[row_ix].split(',')[0].split('|')
        gene_locus = info[2]
        protein_name = info[5]
        gene_name = info[6]
        domain_locus = info[7]
        tooltip_description = '%s[br]Gene Location: %s[br]%s[br]Domain Location: %sAA' %(gene_name, gene_locus, protein_name, domain_locus)
        node_tooltip_description.append(tooltip_description)
    edge_tooltip_description = []
    for x,y in edge_list:
        pair_i_x = node_dict[x]
        pair_i_y = node_dict[y]
        row_index = domain_index[pair_i_x][0]
        col_index = domain_index[pair_i_y][0]
        if row_index < col_index:
            row = DM_input[col_index]
            col = row_index
        else:
            row = DM_input[row_index]
            col = col_index
        score_list = row.split(",")
        score = str(round(1-float(score_list[col + 1]), 3))
        description = 'Domain sequence similarity score by Mafft: [br]%s' %(score)
        edge_tooltip_description.append(description)
    return node_tooltip_description, edge_tooltip_description

def plotting_network(plt_width, nodes, labels, positions, xlim, ylim, edges, legend_x, legend_y, outfile, speciesName_file_list,node_index_dict, KS_color, node_edge_color, nodes_legend, nodes_legend_labels, nodes_legend_pos, nodes_legend_color, nodes_legend_annotation,domain_same_gene, edge_sign_axisX, edge_sign_axisY, speciesName_new_bgc, new_bgc_id):
    labels1 = labels
    labels2 = nodes_legend_labels
    nodes_sum = nodes + nodes_legend
    #print len(nodes_sum), len(nodes), len(nodes_legend)
    G = nx.Graph()
    G.add_nodes_from(nodes_sum)
    edges = edges
    G.add_edges_from(edges)
    fixed_positions = dict(positions, **nodes_legend_pos)
    fixed_nodes = fixed_positions.keys()
    pos = nx.spring_layout(G, pos = fixed_positions, fixed = fixed_nodes)
    #-- plot the defined network
    fig = plt.figure()
    fig.set_size_inches(w = plt_width, h = max(legend_y.values())*1.5)
    ax = fig.add_subplot(111)
    axes = plt.gca()
    axes.set_xlim(xlim)
    axes.set_ylim(ylim)
    #-- map the species name to a list of assembly lines
    species_ID_dict = GSN.get_speciesName(speciesName_file_list = speciesName_file_list, gbID_list = legend_x.keys(), speciesName_new_bgc = speciesName_new_bgc, new_bgc_id = new_bgc_id)
    for i in legend_x.keys():
        ax.annotate(species_ID_dict[i], xy = (legend_x[i], legend_y[i]), xytext = (legend_x[i]-0.2, legend_y[i]), horizontalalignment='right', size = 14)
    for ID in node_index_dict.keys():
        NODES = nx.draw_networkx_nodes(G, pos, nodelist = node_index_dict[ID], node_color = KS_color[ID], node_size = 1000, alpha = 0.7)
        NODES.set_edgecolor(node_edge_color[ID])
    nx.draw_networkx_edges(G, pos, edgelist = edges, width = 8, alpha =0.5, edge_color = '#585858')  # this is to add a line between domains which is a hit in the alignment
    nx.draw_networkx_labels(G, pos, labels = labels1, font_size=11)
    nx.draw_networkx_edges(G, pos, edgelist = domain_same_gene, width = 2, alpha = 0.7, edge_color ='#cd5c5c') # this is to add a hyphen between domains on the same gene
    for l in nodes_legend:
        ax.annotate(nodes_legend_annotation[l], xy = nodes_legend_pos[l], xytext = (nodes_legend_pos[l][0]+0.2, nodes_legend_pos[l][1]), size = 11, horizontalalignment = 'left')
    nodes_in_legend = nx.draw_networkx_nodes(G, pos, nodelist = nodes_legend, node_color = [nodes_legend_color[n] for n in nodes_legend], node_size = 750, alpha = 0.7)
    nodes_in_legend.set_edgecolor([nodes_legend_color[n] for n in nodes_legend])
    nx.draw_networkx_labels(G, pos, labels = labels2, font_size= 11)
    #-- add a marker of pentagram to indicate the break in the original break
    plt.scatter(edge_sign_axisX, edge_sign_axisY, marker = '*', color = 'r', s = 40)
    plt.axis('off') #do not plot the axis
    plt.savefig(outfile, format = 'svg')
    return outfile

def run_plot_network_BGC_list(dist_df, new_bgc_id, new_cluster, new_cluster_index, annotateKS_training,
        annotateKS_novel, speciesName_file_list, dist, color_code_file, outfile, plt_width,
        KS_info_file_list, genomeSize_file_list, ksa_per_new_cluster, genomeSize_new_bgc,
        speciesName_new_bgc, seq_simScore, svg_fout):
    logging.info("TRANSATPKS: creating figure of domain alignment")
    BGC_list_orig = get_bgc_dendrogram_upgma(dist_df)
    global line_gap, feature_loc, feature_edge, feature_strand, edge_len
    feature_loc = 6
    feature_edge = 2
    feature_strand = 3
    line_gap = 0.7
    edge_len = 10000
    novel = []
    training = []

    annoKS_t = [t for t in open(annotateKS_training, 'r').read().split('\n') if t != '']
    annoKS_n1 = [t for t in open(annotateKS_novel, 'r').read().split('\n') if t != '']
    annoKS_n2 = [re.sub('clade_not_conserved', 'noInfo', n) for n in annoKS_n1]
    annoKS_n = [re.sub('out_group_like', 'noInfo', n) for n in annoKS_n2]
    annoKS = annoKS_t[1:] + annoKS_n[1:]
    KSnames = {}
    for a in annoKS:
        KSnames[a.split('\t')[0]] = a
    
    ## The below is added to guard against incompletely annotated pathways which are not included in annoKS.
    ## if statement will fail when the pathway in question is either incomplete or the query sequence
    ## Note: this can subset down the # of pathways to less than the cutoff flag, so if there are less,
    ## this is likely why
    ## Note2: Updated implementation has incomplete pathways included now
    BGC_list = []
    for b in BGC_list_orig:
        if b in KSnames.keys() or b == new_bgc_id[0]:
            BGC_list.append(b)
    for i in BGC_list:
        if len(i.split('|')[1]) > 4: # This is inherited from previous versions, need to doublecheck relevance
            training.append(i)
        else:
            novel.append(i)

    KS_anno_novel, KS_index_novel = extract_KS_anno_unique(id_list=novel, new_bgc_id=new_bgc_id,
            new_cluster=new_cluster, new_cluster_index=new_cluster_index, annoKS=annoKS,
            annoKS_n=annoKS_n, KSnames=KSnames)
    KS_anno_training, KS_index_training = extract_KS_anno_unique(id_list=training, new_bgc_id=new_bgc_id,
            new_cluster=new_cluster, new_cluster_index=new_cluster_index, annoKS=annoKS,
            annoKS_n=annoKS_n, KSnames=KSnames)

    KS_anno = dict(KS_anno_training, **KS_anno_novel)
    KS_index = dict(KS_index_training, **KS_index_novel)
    seq_ordered = BGC_list
    nodes, labels, positions, xlim, ylim, edge_output, legend_pos_x, legend_pos_y, pos_node, node_index_dict, node_anno_dict = get_node_edge_multiAlignment(seq_ordered = seq_ordered, KS_anno = KS_anno, KS_index = KS_index, dist = dist)
    color_dict, node_edge_color_dict = color_nodes_phyloClade(KSanno = KS_anno, color_code_file = color_code_file)
    nodes_legend, nodes_legend_labels, nodes_legend_pos, nodes_legend_color, nodes_legend_annotation = get_legend(KSanno = KS_anno, color_code_file = color_code_file, pos_node = pos_node, nodes_max = max(nodes)+1)
    node_tooltip_description, link_tooltip_description = get_tooltip_links_nodes(node_dict=node_anno_dict, edge_list=edge_output, seq_simScore=seq_simScore)
    domain_loc_edge = MCI.run_get_moduleLoc_contigEdge(KSindex = KS_index,
            KS_info_file = KS_info_file_list, genomeSize_file=genomeSize_file_list,
            feature_loc = feature_loc, feature_edge = feature_edge,
            feature_strand = feature_strand, edge_len = edge_len,
            ksa_per_new_cluster = ksa_per_new_cluster, genomeSize_new_bgc = genomeSize_new_bgc,
            new_bgc_id = new_bgc_id)
    domain_same_gene, edge_sign_axisX, edge_sign_axisY = get_figureInfo_domain_loc_edge(node_index_dict = node_index_dict,
            positions = positions, domain_loc_edge = domain_loc_edge, KSindex = KS_index)
    svg_file = plotting_network(plt_width = plt_width, nodes = nodes, labels = labels, positions = positions,
            xlim = xlim, ylim = ylim, edges = edge_output, legend_x = legend_pos_x,
            legend_y = legend_pos_y, outfile = outfile, speciesName_file_list = speciesName_file_list,
            node_index_dict = node_index_dict, KS_color = color_dict, node_edge_color = node_edge_color_dict,
            nodes_legend = nodes_legend, nodes_legend_labels = nodes_legend_labels,
            nodes_legend_pos = nodes_legend_pos, nodes_legend_color = nodes_legend_color,
            nodes_legend_annotation = nodes_legend_annotation, domain_same_gene = domain_same_gene,
            edge_sign_axisX = edge_sign_axisX, edge_sign_axisY = edge_sign_axisY,
            speciesName_new_bgc = speciesName_new_bgc, new_bgc_id=new_bgc_id)
    modify_svg(f=svg_file, node_tooltip_description=node_tooltip_description, link_tooltip_description=link_tooltip_description , f_out=svg_fout, cluster_id=re.sub('c','', new_bgc_id[0].split('|')[1]))
