import logging
from Bio.SeqFeature import FeatureLocation, CompoundLocation, SeqFeature
from BCBio import GFF


def check_gff_suitability(options, sequences):
    if not options.gff3:
        return

    try:
        examiner = GFF.GFFExaminer()
        gff_data = examiner.available_limits(open(options.gff3))

        # Check if at least one GFF locus appears in sequence
        gff_ids = set([n[0] for n in gff_data['gff_id']])

        if len(gff_ids) == 1 and len(options.all_record_ids) == 1:
            # If both inputs only have one record, assume is the same, but first check coordinate compatibility
            logging.info("GFF3 and sequence have only one record. Assuming is "
                         "the same as long as coordinates are compatible.")
            limit_info = dict(gff_type=['CDS'])

            record_iter = GFF.parse(open(options.gff3), limit_info=limit_info)
            record = next(record_iter)

            coord_max = max([n.location.end.real for n in record.features])
            if coord_max > len(sequences[0]):
                logging.error('GFF3 record and sequence coordinates are not compatible.')
                raise ValueError('Incompatible GFF record and sequence coordinates')
            else:
                options.single_entries = True

        elif len(gff_ids.intersection(options.all_record_ids)) == 0:
            logging.error('No GFF3 record IDs match any sequence record IDs.')
            raise ValueError("GFF3 record IDs don't match sequence file record IDs.")

        else:
            options.single_entries = False

        # Check GFF contains CDSs
        if not ('CDS',) in gff_data['gff_type']:
            logging.error('GFF3 does not contain any CDS.')
            raise ValueError("No CDS features in GFF3 file.")

        # Check CDS are childless but not parentless
        if 'CDS' in set([n for key in examiner.parent_child_map(open(options.gff3)) for n in key]):
            logging.error('GFF3 structure is not suitable. CDS features must be childless but not parentless.')
            raise ValueError('GFF3 structure is not suitable.')

    except AssertionError as e:
        logging.error('Parsing %r failed: %s', options.gff3, e)
        raise


def run(sequence, options):
    handle = open(options.gff3)
    # If there's only one sequence in both, read all, otherwise, read only appropriate part of GFF3.
    if options.single_entries:
        limit_info = False
    else:
        limit_info = dict(gff_id=[sequence.id])

    for record in GFF.parse(handle, limit_info=limit_info):
        for feature in record.features:
            if feature.type == 'CDS':
                new_features = [feature]
            else:
                new_features = check_sub(feature, sequence)
                if not new_features:
                    continue

            name = feature.id
            if len(name) > 40:
                raise ValueError("Feature ID too long, < 40 characters required: %s" % name)
            locus_tag = None
            if "locus_tag" in feature.qualifiers:
                locus_tag = feature.qualifiers["locus_tag"]

            for qtype in ["gene", "name", "Name"]:
                if qtype in feature.qualifiers:
                    name_tmp = feature.qualifiers[qtype][0]
                    # Assume name/Name to be sane if they don't contain a space
                    if " " in name_tmp:
                        continue
                    name = name_tmp
                    break

            for i, n in enumerate(new_features):
                variant = name
                if len(new_features) > 1:
                    variant = "{0}_{1}".format(name, i)
                n.qualifiers['gene'] = [variant]
                if locus_tag is not None:
                    n.qualifiers["locus_tag"] = locus_tag
                sequence.features.append(n)


def check_sub(feature, sequence):
    new_features = []
    loc_list = []
    qual_list = {}
    topop = []
    for sub in feature.sub_features:
        if sub.sub_features:  # If there are sub_features, go deeper
            new_features.extend(check_sub(sub, sequence))
        elif sub.type == 'CDS':
            loc = [sub.location.start.real, sub.location.end.real]
            loc_list.append(FeatureLocation(loc[0], loc[1], strand=sub.strand))
            # For split features (CDSs), the final feature will have the same qualifiers as the children ONLY if
            # they're the same, i.e.: all children have the same "protein_ID" (key and value).
            for qual in sub.qualifiers.keys():
                if qual not in qual_list:
                    qual_list[qual] = sub.qualifiers[qual]
                if qual in qual_list and not qual_list[qual] == sub.qualifiers[qual]:
                    topop.append(qual)

    for n in topop:  # Pop mismatching qualifers over split features
        qual_list.pop(n, None)
    qual_list.pop('Parent', None)  # Pop parent.

    # Only works in tip of the tree, when there's no new_feature built yet. If there is,
    # it means the script just came out of a check_sub and it's ready to return.
    if not new_features:
        if len(loc_list) > 1:
            loc_list = sorted(loc_list, key=lambda x: x.start.real)
            if loc_list[0].strand == 1:
                new_loc = CompoundLocation(loc_list)
            else:
                new_loc = CompoundLocation(list(reversed(loc_list)))
        elif len(loc_list) == 0:
            return new_features
        else:
            new_loc = loc_list[0]

        new_feature = SeqFeature(new_loc)
        new_feature.qualifiers = qual_list
        new_feature.type = 'CDS'
        trans = new_feature.extract(sequence.seq).translate(stop_symbol='')
        new_feature.qualifiers['translation'] = [str(trans)]
        new_features.append(new_feature)

    return new_features
