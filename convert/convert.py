#!/usr/bin/python3

import re
import sys
from collections import OrderedDict
from getopt import getopt
from glob import glob
from os.path import basename
from os.path import splitext

MALE = '1'
FEMALE = '2'

case_dir = 'cases'
control_dir = 'controls'
unknown_dir = 'unknowns'
recursive = False
family_id = 'FAM001'
spacing = 0
out_dir = '.'
want_parents = True
want_sexes = True

optlist, args = getopt(
    sys.argv[1:], '-ar',
    ['cases=', 'controls=', 'unknowns=', 'recursive', 'family=', 'no-parents', 'no-sexes', 'spacing=', 'out=', 'help']
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
    elif name == '--no-parents':
        want_parents = False
    elif name == '--no-sexes':
        want_sexes = False
        want_parents = False
    elif name == '--out':
        out_dir = value
    elif name == '--spacing':
        spacing = float(value)
    elif name == '--help':
        print('Syntax: ./convert.py [--cases=<dir>] [--controls=<dir>] [--unknowns=<dir>]')
        print('                [-r | --recursive] [--family=<name>] [--no-parents]')
        print('                [--no-sexes] [--spacing=<cm>] [--out=<dir>]')
        print('cases defaults to ./cases, controls defaults to ./controls, unknowns defaults')
        print('to ./unknowns, family defaults to FAM001, and out defaults to the current')
        print('directory.')
        exit()

print('Loading HapMap...')

hap_map = {}

for chromosome in list(map(str, range(1, 23))) + ['X']:
    try:
        f = open('hapmap/genetic_map_GRCh37_chr' + chromosome + '.txt')
        hap_map[chromosome] = [(0, 0)]
        for line in f:
            if line.startswith('Chromosome'):
                continue

            cells = line.split('\t')
            bp_pos = int(cells[1])
            cm_pos = float(cells[3])

            hap_map[chromosome].append((bp_pos, cm_pos))
    except FileNotFoundError:
        continue

print('Loading 23andMe raw data files...')

raw_data = OrderedDict()
affection_table = {}
snp_map = OrderedDict()

def load_files(txt_dir, affection):
    for filename in sorted(glob(txt_dir + '/**', recursive=recursive)):
        person_id = re.sub('\W', '', splitext(basename(filename))[0])
        affection_table[person_id] = affection

        for line in open(filename, 'r'):
            if line.startswith('#'):
                continue

            cells = line.split('\t')
            rsid = cells[0]
            chromosome = cells[1]
            bp_pos = int(cells[2])
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

            if not chromosome in hap_map:
                cm_pos = 0
            else:
                low = 0
                high = len(hap_map[chromosome]) - 1
                while high - low > 1:
                    middle = int((low + high) / 2)
                    if hap_map[chromosome][middle][0] > bp_pos:
                        high = middle
                    else:
                        low = middle
                cm_pos = (
                    hap_map[chromosome][low][1] + (hap_map[chromosome][high][1] - hap_map[chromosome][low][1]) *
                    (bp_pos - hap_map[chromosome][low][0]) / (hap_map[chromosome][high][0] - hap_map[chromosome][low][0])
                )

            if not rsid in snp_map:
                snp_map[rsid] = (chromosome, bp_pos, cm_pos)

load_files(case_dir, '2')
load_files(control_dir, '1')
load_files(unknown_dir, '0')

sex_table = {}

if want_sexes:
    print('Inferring sexes...')

    for person_id in raw_data:
        male_snps = 0
        for rsid in raw_data[person_id]['Y']:
            if raw_data[person_id]['Y'][rsid] != ('0', '0'):
                male_snps += 1
        if male_snps / len(raw_data[person_id]['Y']) > 0.1:
            sex_table[person_id] = MALE
        else:
            sex_table[person_id] = FEMALE

parents_table = {}

if want_parents:
    print('Inferring relationships...')

    PARENT_OR_CHILD = 'Parent/Child'
    SIBLING = 'Sibling'

    relationship_table = {}
    for person_id in raw_data:
        relationship_table[person_id] = {}
    for person1_id in raw_data:
        for person2_id in raw_data:
            if person1_id == person2_id or person1_id in relationship_table[person2_id]:
                continue

            half_matches = 0
            full_matches = 0
            possible_matches = 0

            for chromosome in map(str, range(1, 23)):
                if chromosome in hap_map:
                    last_cm_pos = 0
                    sorted_snps = sorted(raw_data[person1_id][chromosome], key=lambda rsid: snp_map[rsid][2])
                    for rsid in sorted_snps:
                        person1_base1 = raw_data[person1_id][chromosome][rsid][0]
                        person1_base2 = raw_data[person1_id][chromosome][rsid][1]
                        person2_base1 = raw_data[person2_id][chromosome][rsid][0]
                        person2_base2 = raw_data[person2_id][chromosome][rsid][1]
                        cm_pos = snp_map[rsid][2]
                        dist_to_last = cm_pos - last_cm_pos
                        if person1_base1 == person2_base1:
                            half_matches += dist_to_last
                            if person1_base2 == person2_base2:
                                full_matches += dist_to_last
                        elif person1_base1 == person2_base2:
                            half_matches += dist_to_last
                            if person1_base2 == person2_base1:
                                full_matches += dist_to_last
                        elif person1_base2 == person2_base1:
                            half_matches += dist_to_last
                            if person1_base1 == person2_base2:
                                full_matches += dist_to_last
                        elif person1_base2 == person2_base2:
                            half_matches += dist_to_last
                            if person1_base1 == person2_base1:
                                full_matches += dist_to_last
                        last_cm_pos = cm_pos
                    possible_matches += last_cm_pos
                else:
                    snps = raw_data[person1_id][chromosome]
                    for rsid in snps:
                        person1_base1 = raw_data[person1_id][chromosome][rsid][0]
                        person1_base2 = raw_data[person1_id][chromosome][rsid][1]
                        person2_base1 = raw_data[person2_id][chromosome][rsid][0]
                        person2_base2 = raw_data[person2_id][chromosome][rsid][1]
                        if person1_base1 == person2_base1:
                            half_matches += 1
                            if person1_base2 == person2_base2:
                                full_matches += 1
                        elif person1_base1 == person2_base2:
                            half_matches += 1
                            if person1_base2 == person2_base1:
                                full_matches += 1
                        elif person1_base2 == person2_base1:
                            half_matches += 1
                            if person1_base1 == person2_base2:
                                full_matches += 1
                        elif person1_base2 == person2_base2:
                            half_matches += 1
                            if person1_base1 == person2_base1:
                                full_matches += 1
                    possible_matches += len(snps)

            '''
            print(person1_id + ' to ' + person2_id + ':' +
                ' half ' + str(half_matches / possible_matches) +
                ' full ' + str(full_matches / possible_matches))
            #'''
            if half_matches / possible_matches > 0.98:
                relationship_table[person1_id][person2_id] = PARENT_OR_CHILD
                relationship_table[person2_id][person1_id] = PARENT_OR_CHILD
            elif full_matches / possible_matches > 0.70:
                relationship_table[person1_id][person2_id] = SIBLING
                relationship_table[person2_id][person1_id] = SIBLING

    #print(relationship_table)

    for proband_id in list(raw_data.keys()):
        father_id = '0'
        mother_id = '0'
        for potential_parent_id in relationship_table[proband_id]:
            if relationship_table[proband_id][potential_parent_id] != PARENT_OR_CHILD:
                continue
            for sibling_id in relationship_table[proband_id]:
                if relationship_table[proband_id][sibling_id] != SIBLING:
                    continue
                if (potential_parent_id in relationship_table[sibling_id] and
                        relationship_table[sibling_id][potential_parent_id] == PARENT_OR_CHILD):
                    if sex_table[potential_parent_id] == MALE:
                        father_id = potential_parent_id
                    else:
                        mother_id = potential_parent_id
                    if mother_id != '0' and father_id != '0':
                        break
            if mother_id != '0' and father_id != '0':
                break

        if father_id == '0' and mother_id != '0':
            missing_parent_id = mother_id + 'Spouse'
            missing_parent_sex = MALE
            father_id = missing_parent_id
        elif mother_id == '0' and father_id != '0':
            missing_parent_id = father_id + 'Spouse'
            missing_parent_sex = FEMALE
            mother_id = missing_parent_id
        else:
            missing_parent_id = None

        parents_table[proband_id] = (father_id, mother_id)

        if missing_parent_id and missing_parent_id not in raw_data:
            affection_table[missing_parent_id] = '0'
            raw_data[missing_parent_id] = {}
            for rsid in snp_map:
                chromosome = snp_map[rsid][0]
                if chromosome not in raw_data[missing_parent_id]:
                    raw_data[missing_parent_id][chromosome] = {}
                raw_data[missing_parent_id][chromosome][rsid] = ('0', '0')
            sex_table[missing_parent_id] = missing_parent_sex
            parents_table[missing_parent_id] = ('0', '0')

if spacing > 0:
    print('Thinning the data...')

    new_snp_map = OrderedDict()
    prev_chromosome = None
    next_cm_pos = 0

    for rsid in snp_map:
        if snp_map[rsid][0] != prev_chromosome:
            next_cm_pos = 0
            prev_chromosome = snp_map[rsid][0]

        if snp_map[rsid][2] >= next_cm_pos:
            new_snp_map[rsid] = snp_map[rsid]
            next_cm_pos += spacing

    snp_map = new_snp_map

print('Writing files...')

ped_file = open(out_dir + '/' + family_id + '.ped', 'w')
for proband_id in raw_data:
    if want_parents:
        father_id = parents_table[proband_id][0]
        mother_id = parents_table[proband_id][1]
    else:
        father_id = '0'
        mother_id = '0'

    if want_sexes:
        sex = sex_table[proband_id]
    else:
        sex = '0'

    ped_file.write(
        family_id + '\t' +
        proband_id + '\t' +
        father_id + '\t' +
        mother_id + '\t' +
        sex + '\t' +
        affection_table[proband_id]
    )
    for rsid in snp_map:
        chromosome = snp_map[rsid][0]
        bases = raw_data[proband_id][chromosome][rsid]
        ped_file.write('\t' + bases[0] + ' ' + bases[1])
    ped_file.write('\n')

plink_map_file = open(out_dir + '/' + family_id + '.plink.map', 'w')
plink_map_file.write('# chromosome\trsid\tmorgans\tposition\n')
for rsid in snp_map:
    plink_map_file.write(snp_map[rsid][0] + '\t' + rsid + '\t' + str(snp_map[rsid][2] / 100) + '\t' + str(snp_map[rsid][1]) + '\n')

merlin_map_file = open(out_dir + '/' + family_id + '.merlin.map', 'w')
merlin_map_file.write('CHROMOSOME\tMARKER_NAME\tPOSITION\n')
for rsid in snp_map:
    merlin_map_file.write(snp_map[rsid][0] + '\t' + rsid + '\t' + str(snp_map[rsid][2]) + '\n')

merlin_dat_file = open(out_dir + '/' + family_id + '.merlin.dat', 'w')
merlin_dat_file.write('A  affection\n')
for rsid in snp_map:
    merlin_dat_file.write('M  ' + rsid + '\n')
