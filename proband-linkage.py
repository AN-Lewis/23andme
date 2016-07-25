#!/usr/bin/python3
#
# 23andMe Proband Linkage Analyzer
#
# Copyright (c) 2016 Alex Henrie
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import csv
import glob
import sys
from copy import copy
from getopt import getopt
from math import sqrt
from os.path import basename
from os.path import splitext
from scipy.stats import norm
from sys import stderr

case_dir = 'cases'
control_dir = 'controls'
proband_is_case = None
alpha = 0.05
m = 23
want_misfits = False

optlist, args = getopt(sys.argv[1:], '-a', ['cases=', 'controls=', 'proband=', 'alpha=', 'no-bonferroni', 'misfits', 'help'])
for name, value in optlist:
    if name == '--cases':
        case_dir = value
    elif name == '--controls':
        control_dir = value
    elif name == '--proband':
        if value.lower() == 'case':
            proband_is_case = True
        elif value.lower() == 'control':
            proband_is_case = False
    elif name in ('-a', '--alpha'):
        alpha = float(value)
    elif name == '--no-bonferroni':
        m = 1
    elif name == '--misfits':
        want_misfits = True
    elif name == '--help':
        print('Syntax: ./proband-linkage.py [--cases=<dir>] [--controls=<dir>] [--proband=(case|control)] [--alpha=<value>] [--no-bonferroni] [--misfits]')
        print('cases defaults to ./cases, controls defaults to ./controls, and alpha defaults to 0.05.')
        exit(1)

while proband_is_case == None:
    stderr.write('Is the proband a case or a control? ')
    proband_is_case = input()
    if proband_is_case.lower() == 'case':
        proband_is_case = True
    elif proband_is_case.lower() == 'control':
        proband_is_case = False
    else:
        proband_is_case = None

def create_start_point(segments, point):
    for i in range(0, len(segments)):
        if segments[i][0] < point and point < segments[i][1]:
            #shrink the segment and insert a new segment after it with the same people
            segments.insert(i + 1, [point, segments[i][1], copy(segments[i][2]), copy(segments[i][3])])
            segments[i][1] = point - 1
            break

def create_end_point(segments, point):
    for i in range(0, len(segments)):
        if segments[i][0] < point and point < segments[i][1]:
            #shrink the segment and insert a new segment before it with the same people
            segments.insert(i, [segments[i][0], point, copy(segments[i][2]), copy(segments[i][3])])
            segments[i + 1][0] = point + 1
            break

if proband_is_case:
    cases = ['Proband']
    controls = []
else:
    cases = []
    controls = ['Proband']

linkage = {
     '1': [[0, sys.maxsize, copy(cases), copy(controls)]],
     '2': [[0, sys.maxsize, copy(cases), copy(controls)]],
     '3': [[0, sys.maxsize, copy(cases), copy(controls)]],
     '4': [[0, sys.maxsize, copy(cases), copy(controls)]],
     '5': [[0, sys.maxsize, copy(cases), copy(controls)]],
     '6': [[0, sys.maxsize, copy(cases), copy(controls)]],
     '7': [[0, sys.maxsize, copy(cases), copy(controls)]],
     '8': [[0, sys.maxsize, copy(cases), copy(controls)]],
     '9': [[0, sys.maxsize, copy(cases), copy(controls)]],
    '10': [[0, sys.maxsize, copy(cases), copy(controls)]],
    '11': [[0, sys.maxsize, copy(cases), copy(controls)]],
    '12': [[0, sys.maxsize, copy(cases), copy(controls)]],
    '13': [[0, sys.maxsize, copy(cases), copy(controls)]],
    '14': [[0, sys.maxsize, copy(cases), copy(controls)]],
    '15': [[0, sys.maxsize, copy(cases), copy(controls)]],
    '16': [[0, sys.maxsize, copy(cases), copy(controls)]],
    '17': [[0, sys.maxsize, copy(cases), copy(controls)]],
    '18': [[0, sys.maxsize, copy(cases), copy(controls)]],
    '19': [[0, sys.maxsize, copy(cases), copy(controls)]],
    '20': [[0, sys.maxsize, copy(cases), copy(controls)]],
    '21': [[0, sys.maxsize, copy(cases), copy(controls)]],
    '22': [[0, sys.maxsize, copy(cases), copy(controls)]],
     'X': [[0, sys.maxsize, copy(cases), copy(controls)]],
}

def load_files(csv_dir, is_case):
    global cases
    global controls
    for filename in glob.glob(csv_dir + '/*.csv'):
        person_name = splitext(basename(filename))[0]
        for row in csv.reader(open(filename, 'r')):
            if row == ['Comparison', 'Chromosome', 'Start Point', 'End Point', 'Genetic Distance', '#SNPs']:
                continue
            chromosome = row[1]
            start_point = int(row[2])
            end_point = int(row[3])
            create_start_point(linkage[chromosome], start_point)
            create_end_point(linkage[chromosome], end_point)
            for segment in linkage[chromosome]:
                if start_point <= segment[0] and segment[1] <= end_point:
                    if is_case:
                        segment[2].append(person_name)
                    else:
                        segment[3].append(person_name)
        if is_case:
            cases.append(person_name)
        else:
            controls.append(person_name)

load_files(case_dir, True)
load_files(control_dir, False)

total_segments = 0
for chromosome, segments in linkage.items():
    total_segments += len(segments)

stderr.write('Comparing ' + str(total_segments) + ' segments between ' + str(len(cases)) + ' cases and ' + str(len(controls)) + ' controls.\n')

if len(cases) < 10 or len(controls) < 10:
    stderr.write('WARNING: You have less than 10 cases or less than 10 controls.\n')
    stderr.write('The p values output by this program are not accurate.\n')

difference_found = False

def pad_key(pair):
    try:
        return (str(int(pair[0])).zfill(2), pair[1])
    except ValueError:
        return pair

for chromosome, segments in sorted(linkage.items(), key=pad_key):
    for segment in segments:
        if segment[0] == 0 or segment[1] == sys.maxsize:
            continue
        try:
            p = (len(segment[2]) + len(segment[3])) / (len(cases) + len(controls))
            se = sqrt(p * (1 - p) * (1 / len(cases) + 1 / len(controls)))
            z = -abs((len(segment[2]) / len(cases) - len(segment[3]) / len(controls)) / se)
            p = 2 * norm.cdf(z)
            if p <= alpha / m:
                if not difference_found:
                    header = 'Chromosome\tStart\tEnd\tCase freq\tControl freq\tp'
                    if want_misfits:
                        header += '\tCase misfits\tControl misfits'
                    print(header)
                    difference_found = True

                case_freq = len(segment[2]) / len(cases)
                control_freq = len(segment[3]) / len(controls)
                if case_freq >= control_freq:
                    #we expect the cases to have the allele and the controls to not have it
                    case_misfits = set(cases).difference(set(segment[2]))
                    control_misfits = set(segment[3])
                else:
                    #we expect the cases to not have the allele and the controls to have it
                    case_misfits = set(segment[2])
                    control_misfits = set(controls).difference(set(segment[3]))

                output = (
                    str(chromosome) + '\t' +
                    str(segment[0]) + '\t' +
                    str(segment[1]) + '\t' +
                    str(case_freq) + '\t' +
                    str(control_freq) + '\t' +
                    str(p)
                )
                if want_misfits:
                    output += (
                        '\t' + str(case_misfits) +
                        '\t' + str(control_misfits)
                    )
                print(output)
        except ZeroDivisionError:
            pass

if not difference_found:
    stderr.write('There is no significant difference between the cases and the controls.\n')
