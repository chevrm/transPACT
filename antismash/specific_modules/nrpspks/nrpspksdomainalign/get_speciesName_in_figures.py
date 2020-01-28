#----------------------------------------------------------------------------------------------------------
## Copyright (c) 2016 Xiaowen Lu
## Wageningen University
## Bioinformatics Group
#
#  This is the function to get the species name based on the keys of the output dict from the function
#  extract_KS_anno_unique in '/Users/Xiaowen/Documents/trans_AT/scripts/remove_same_transATPKS.py'
#----------------------------------------------------------------------------------------------------------


# speciesName_file_list = ['/Users/Xiaowen/Documents/trans_AT/KS_domains_PredtransAT_output_20160330/species_ID_training_transATPKS.txt', '/Users/Xiaowen/Documents/trans_AT/KS_domains_PredtransAT_output_20160330/species_ID_pred_transATPKS.txt']
# gbID_list = ['CGFR01000002|c1', 'AJST01000001|c6(extra 2)', 'CXWA01000005|c1']


def get_speciesName(speciesName_file_list, gbID_list, speciesName_new_bgc, new_bgc_id):

    speciesName = {}
    for l in speciesName_file_list:
        infile = [s for s in open(l, 'r').read().split('\n') if s != '']
        for i in infile:
            k = i.split('\t')[0]
            v = i.split('\t')[1]
            speciesName[k] = v

    for bgc in new_bgc_id:
        speciesName[bgc] = '|'.join((speciesName_new_bgc, bgc))


    species_output = {}
    for id in gbID_list:
        id_m = id.split('(')
        species = speciesName[id_m[0]]
        if len(id_m) == 1:
            species_output[id] = species
        else:
            species_output[id] = '('.join([species, id_m[1]])

    return species_output

