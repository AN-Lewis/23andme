#!/usr/bin/python3

import csv
from collections import OrderedDict
from sys import argv

rows_to_ignore = [
    ['rsid', 'chromosome', 'position', 'allele1', 'allele2'], #Ancestry.com
    ['Name', 'Variation', 'Chromosome', 'Position', 'Strand', 'YourCode'], #deCODEme
    ['RSID', 'CHROMOSOME', 'POSITION', 'RESULT'], #FTDNA
]

all_snps = set()
snps_by_file = []

for filename in argv[1:]:
    snps_in_file = set()
    if filename.endswith('.csv'):
        delimiter = ','
    else:
        delimiter = '\t'
    for row in csv.reader(open(filename, 'r'), delimiter=delimiter):
        if row[0].startswith('#'):
            continue
        if row in rows_to_ignore:
            continue
        rsid = row[0]
        all_snps.add(rsid)
        snps_in_file.add(rsid)
    snps_by_file.append(snps_in_file)

print('\t' + '\t'.join(argv[1:]))
for rsid in sorted(all_snps):
    availability = []
    for snps_in_file in snps_by_file:
        if rsid in snps_in_file:
            availability.append('Yes')
        else:
            availability.append('No')
    print(rsid + '\t' + '\t'.join(availability))
