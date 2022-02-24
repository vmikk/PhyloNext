#!/usr/bin/python
"""
Example command (in python3):
python fetch_gbif_extinct_tax.py 359 --outfile extinct_tax_ids.txt

highertaxon_key should be GBIF taxon key, e.g. 358 for Reptilia, 359 for Mammalia, etc.
"""

import argparse
import json
import requests
import sys

parser = argparse.ArgumentParser(description="Python script to fetch extinct taxon ids from GBIF API")
parser.add_argument("highertaxon_key", help="Need highertaxon_key in numeric format!")
parser.add_argument("--outfile", default="extinct_tax_ids.txt")
args = parser.parse_args()

# read in args
highertaxon_key = args.highertaxon_key
if not highertaxon_key.isdigit():
    raise ValueError("highertaxon_key has to be numeric", highertaxon_key)

outfile = args.outfile

# Search parameters (target GBIF backbone taxonomy)
rank = "SPECIES"
status = "ACCEPTED"
is_extinct = True
dataset_key = "d7dddbf4-2cf0-4f39-9b2a-bb099caae36c"
limit = 10
offset = 0

# Fetch the data
base_query = "https://api.gbif.org/v1/species/search"

# Get record count
response_ct = requests.get(base_query, params={
    "rank": f"{rank}",
    "highertaxon_key": f"{highertaxon_key}",
    "status": f"{status}",
    "datasetKey": f"{dataset_key}",
    "isExtinct": f"{is_extinct}",
    "limit": "0",
    "offset": "0",
})
if response_ct.status_code == 200:
    contents_ct = response_ct.json()

    if contents_ct["count"] > 0:
        print("No. of records to be fetched: " + str(contents_ct["count"]))

        record_ct = 0
        proceed = True
        
        # Start fetching the data
        with open(outfile, "w") as o:
            while proceed:
                print(f"Starting to add {limit} to already downloaded {record_ct}")
                response = requests.get(base_query, params={
                    "rank": f"{rank}",
                    "highertaxon_key": f"{highertaxon_key}",
                    "status": f"{status}",
                    "datasetKey": f"{dataset_key}",
                    "isExtinct": f"{is_extinct}",
                    "limit": f"{limit}",
                    "offset": record_ct,
                })
                if response.status_code == 200:
                    contents = response.json()
                    results = contents["results"]
                    tax_key_list = [str(record["key"]) for record in results]
                    tax_key_str = ",".join(tax_key_list)
                    o.write(tax_key_str)

                    total_count = len(results)
                    record_ct += total_count
                    if total_count < limit:
                        proceed = False
                    else:
                        o.write(",")
                else:
                    sys.exit("Command exited with status code " + str(response.status_code))
