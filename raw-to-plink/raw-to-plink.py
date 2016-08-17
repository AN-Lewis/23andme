#!/usr/bin/python3

import re
import sys
from collections import OrderedDict
from getopt import getopt
from glob import glob
from os.path import basename
from os.path import splitext

case_dir = 'cases'
control_dir = 'controls'
unknown_dir = 'unknowns'
family_id = 'FAM001'
recursive = False

optlist, args = getopt(
    sys.argv[1:], '-ar',
    ['cases=', 'controls=', 'unknowns=', 'recursive', 'family', 'help']
)
for name, value in optlist:
    if name == '--cases':
        case_dir = value
    elif name == '--controls':
        control_dir = value
    elif name == '--unknowns':
        unknown_dir = value
    elif name in ('-r', '--recursive'):
        recursive = True
    elif name == '--family':
        family_id = value;
    elif name == '--help':
        print('Syntax: ./raw-to-plink.py [--cases=<dir>] [--controls=<dir>] [--unknowns=<dir>]')
        print('                [--family=<name>] [-r | --recursive]')
        print('cases defaults to ./cases, controls defaults to ./controls, unknowns defaults')
        print('to ./unknowns, and family defaults to FAM001.')
        exit()

raw_data = OrderedDict()
affections = {}
snp_map = OrderedDict()

def load_files(txt_dir, affection):
    for filename in sorted(glob(txt_dir + '/**', recursive=recursive)):
        person_id = re.sub('\W', '', splitext(basename(filename))[0])
        affections[person_id] = affection

        for line in open(filename, 'r'):
            if line.startswith('#'):
                continue

            cells = line.split('\t')
            rsid = cells[0]
            chromosome = cells[1]
            position = cells[2]
            base1 = cells[3][0].replace('-', '0')
            base2 = cells[3][1].replace('-', '0')
            if base2 == '\n': #X, Y, MT
                base2 = base1
            if base1 == '0':
                base2 = '0'
            if base2 == '0':
                base1 = '0'

            if not person_id in raw_data:
                raw_data[person_id] = {}
            if not chromosome in raw_data[person_id]:
                raw_data[person_id][chromosome] = {}
            raw_data[person_id][chromosome][rsid] = (base1, base2)

            if not rsid in snp_map:
                snp_map[rsid] = (chromosome, position)

load_files(case_dir, '2')
load_files(control_dir, '1')
load_files(unknown_dir, '0')

MALE = 1
FEMALE = 2
def get_sex(person_id):
    male_snps = 0
    for rsid in raw_data[person_id]['Y']:
        if raw_data[person_id]['Y'][rsid] != ('0', '0'):
            male_snps += 1
    if male_snps / len(raw_data[person_id]['Y']) > 0.1:
        return MALE
    else:
        return FEMALE

ped_file = open(family_id + '.ped', 'w')
for proband_id in raw_data:
    #TODO: automatically infer parent-child relationships
    ped_file.write(
        family_id + '\t' +
        proband_id + '\t' +
        '0\t' +
        '0\t' +
        str(get_sex(proband_id)) + '\t' +
        affections[proband_id]
    )
    for rsid in snp_map:
        chromosome = snp_map[rsid][0]
        bases = raw_data[proband_id][chromosome][rsid]
        if '\n' in bases[0] or '\n' in bases[1]:
            print(rsid)
        ped_file.write('\t' + bases[0] + ' ' + bases[1])
    ped_file.write('\n')

map_file = open(family_id + '.map', 'w')
map_file.write('# chromosome\trsid\tmorgans\tposition\n')
for rsid in snp_map:
    map_file.write(snp_map[rsid][0] + '\t' + rsid + '\t0\t' + snp_map[rsid][1] + '\n')
