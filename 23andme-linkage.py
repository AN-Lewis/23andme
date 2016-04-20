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
from getopt import getopt
from math import exp
from math import pi
from math import sqrt

case_dir = 'cases'
control_dir = 'controls'
proband_is_case = None
alpha = 0.05

optlist, args = getopt(sys.argv[1:], '-a', ['cases=', 'controls=', 'proband=', 'alpha=', 'help'])
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
    elif name == '--help':
        print('Syntax: ./23andme-linkage.py [--cases=<dir>] [--controls=<dir>] [--proband=(case|control)] [--alpha=<value>]')
        print('cases defaults to ./cases, controls defaults to ./controls, and alpha defaults to 0.05.')
        exit(1)

while proband_is_case == None:
    proband_is_case = input('Is the proband a case or a control? ')
    if proband_is_case.lower() == 'case':
        proband_is_case = True
    elif proband_is_case.lower() == 'control':
        proband_is_case = False
    else:
        proband_is_case = None

def create_start_point(segments, point):
    for i in range(0, len(segments)):
        if segments[i][0] < point and point < segments[i][1]:
            #shrink the segment and insert a new segment after it with the same totals
            segments.insert(i + 1, [point, segments[i][1], segments[i][2], segments[i][3]])
            segments[i][1] = point - 1
            break

def create_end_point(segments, point):
    for i in range(0, len(segments)):
        if segments[i][0] < point and point < segments[i][1]:
            #shrink the segment and insert a new segment before it with the same totals
            segments.insert(i, [segments[i][0], point, segments[i][2], segments[i][3]])
            segments[i + 1][0] = point + 1
            break

linkage = {
     '1': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
     '2': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
     '3': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
     '4': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
     '5': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
     '6': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
     '7': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
     '8': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
     '9': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
    '10': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
    '11': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
    '12': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
    '13': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
    '14': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
    '15': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
    '16': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
    '17': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
    '18': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
    '19': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
    '20': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
    '21': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
    '22': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
     'X': [[0, sys.maxsize, int(proband_is_case), int(not proband_is_case)]],
}

cases = int(proband_is_case)
controls = int(not proband_is_case)

def load_files(csv_dir, is_case):
    global cases
    global controls
    for filename in glob.glob(csv_dir + '/*.csv'):
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
                        segment[2] += 1
                    else:
                        segment[3] += 1
        if is_case:
            cases += 1
        else:
            controls += 1

load_files(case_dir, True)
load_files(control_dir, False)

print('Comparing ' + str(cases) + ' cases to ' + str(controls) + ' controls.')

if cases < 10 or controls < 10:
    print('WARNING: You do not have enough cases and controls.')
    print('The p values output by this program are not accurate.')

difference_found = False

for chromosome, segments in linkage.items():
    for segment in segments:
        if segment[0] == 0 or segment[1] == sys.maxsize:
            continue
        try:
            p = (segment[2] + segment[3]) / (cases + controls)
            se = sqrt(p * (1 - p) * (1 / cases + 1 / controls))
            z = (segment[2] / cases - segment[3] / controls) / se
            p = 1 - exp(-z**2 / 2) / sqrt(2 * pi)
            if p <= alpha:
                if not difference_found:
                    print('Chromosome\tStart\tEnd\tCase freq\tControl freq\tp')
                    difference_found = True
                print(
                    str(chromosome) + '\t' +
                    str(segment[0]) + '\t' +
                    str(segment[1]) + '\t' +
                    str(segment[2] / cases) + '\t' +
                    str(segment[3] / controls) + '\t' +
                    str(p)
                )
        except ZeroDivisionError:
            pass

if not difference_found:
    print('There is no significant difference between the cases and the controls.')
