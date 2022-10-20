
#!/usr/bin/env python3

import argparse
import os
import dendropy

from chronosynth import chronogram
from opentree import OT, annotations, taxonomy_helpers

from opentree.taxonomy_helpers import labelled_induced_synth




parser = argparse.ArgumentParser()
parser.add_argument("-q","--query", help="File containing ott_ids. First should have the header 'ott_id' and contain OpenTree taxon ids")
parser.add_argument("-ma","--max_age", help="Root age estimate")
parser.add_argument("-o","--output", default="bladj.tre", help="output file name")
parser.add_argument("-l","--label_format", default="name_and_id",help="label format. One of 'name', 'id', 'name_and_id'. Default is 'name_and_id'")

args = parser.parse_args()




taxa = {}

assert os.path.exists(args.query)

queryfile = open(args.query)
header = queryfile.readline().strip().split(',')

for lin in queryfile:
    lii = lin.strip().split(',')
    taxa['ott'+lii[0]]=dict(zip(header, lii))


queryfile.close()

ott_ids = list(taxa.keys())



chronogram.date_synth_subtree(node_ids=ott_ids,
                              max_age=args.max_age,
                              output_dir='.',
                              method="bladj",
                              summary=args.output)


#resp = OT.synth_induced_tree(node_ids = ott_ids, label_format = 'id', ignore_unknown_ids=True)


tree = dendropy.Tree.get_from_path(args.output, schema= "newick")

node_annotations = annotations.generate_synth_node_annotation(tree)

annotations.write_itol_conflict(node_annotations)
annotations.write_itol_support(node_annotations)



translation_dict = {}
for node in node_annotations:
    if node in taxa:
        translation_dict[node]=taxa[node]['unique_name']

annotations.write_itol_relabel(translation_dict, filename="ottlabel.txt")
