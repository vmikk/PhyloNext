#!/usr/bin/env python
"""
Script to get an induced subtree from synthetic phylogenetic tree for a set of OpenTree IDs
Author:  Emily Jane McTavish
Version: b1ead89 (September 20, 2022)
Link:    https://github.com/McTavishLab/GBIF-Biodiverse-OpenTree/

to run:
python induced_synth_subtree_from_csv.py --query ../tests/query.csv --output_dir amph_tree

-ott-ids data file should have a column containing ott_ids, labeled 'ott_id', and be comma delimited

Will generate a tree, a logfile and a citations file, all in the output directory.
"""

import argparse
import os
import json
import requests
import dendropy
from opentree.taxonomy_helpers import labelled_induced_synth
from opentree import OT, annotations


parser = argparse.ArgumentParser()
parser.add_argument(
    "-q",
    "--query",
    help="File containing ott_ids. First column should have the header 'ott_id' and contain OpenTree taxon ids",
)
parser.add_argument(
    "-o",
    "--output_dir",
    default="synth_output",
    help="output folder. default is synth_dir",
)
parser.add_argument(
    "-l",
    "--label_format",
    default="name_and_id",
    help="label format. One of 'name', 'id', 'name_and_id'. Default is 'name_and_id'",
)

args = parser.parse_args()

print("Running script 'induced_synth_subtree_from_csv.py'\n", file=sys.stderr)

assert os.path.exists(args.query)

if not os.path.exists(args.output_dir):
    os.mkdir(args.output_dir)

queryfile = open(args.query)
print("Input query:", args.query, "\n", file=sys.stderr)
header = queryfile.readline()
assert header.split(",")[0] == "ott_id"


taxa = {}

assert os.path.exists(args.query)

queryfile = open(args.query)
header = queryfile.readline().strip().split(",")

for lin in queryfile:
    lii = lin.strip().split(",")
    taxa["ott" + lii[0]] = dict(zip(header, lii))


queryfile.close()

node_ids = list(taxa.keys())
ott_ids = [node.strip("ott") for node in node_ids]

print("Preparing synth tree\n", file=sys.stderr)
ret = labelled_induced_synth(
    ott_ids=ott_ids,
    label_format=args.label_format,
    inc_unlabelled_mrca=False,
    standardize=True,
)

print("Exporting labelled tree\n", file=sys.stderr)
ret["labelled_tree"].write(path=args.output_dir + "/labelled_tree.tre", schema="newick")

log = open("{}/synth.log".format(args.output_dir), "w")
log.write("unknown_ids:\n")
log.write("\n".join(ret["unknown_ids"]))
log.write("\nnon-monophyletic_taxa:\n")
log.write("\n".join(ret["non-monophyletic_taxa"]))
log.close()


cites_file = open("{}/citations.txt".format(args.output_dir), "w")
citations = OT.get_citations(ret["supporting_studies"])

cites_file.write(citations)


session = requests.Session()

url = "https://dates.opentreeoflife.org/v4/dates/dated_tree"
payload = {"node_ids": node_ids}

resp = session.post(url=url, data=json.dumps(payload))
resp_dict = json.loads(resp.content.decode())

print("Obtaining dated tree\n", file=sys.stderr)
dated_tree = dendropy.Tree.get(
    data=resp_dict["dated_trees_newick_list"][0], schema="newick"
)

print("Exporting dated tree\n", file=sys.stderr)
dated_tree.write(path=args.output_dir + "/ottid_dated_tree.tre", schema="newick")


date_cites_file = open("{}/date_citations.txt".format(args.output_dir), "w")
date_citations = OT.get_citations(ret["supporting_studies"])

date_cites_file.write(date_citations)

print("Node annotations\n", file=sys.stderr)
node_annotations = annotations.generate_synth_node_annotation(dated_tree)

translation_dict = {}
for node in node_annotations:
    if node in taxa:
        translation_dict[node] = taxa[node]["unique_name"]

annotations.write_itol_relabel(
    translation_dict, filename=args.output_dir + "/ottlabel.txt"
)

for tip in dated_tree.leaf_node_iter():
    tip.taxon.label = translation_dict[tip.taxon.label]


dated_tree.write(path=args.output_dir + "/ott_label_dated_tree.tre", schema="newick")


annotations.write_itol_conflict(
    node_annotations,
    filename=args.output_dir + "/conflict_annot.tre",
)
annotations.write_itol_support(
    node_annotations,
    filename=args.output_dir + "/support_annot.tre",
)

print("All done.\n", file=sys.stderr)
