#!/usr/bin/env python
"""
:File: goshifter.py
:Author: Gosia Trynka <gosia@broadinstitute.org>
:Last updated: 2014-07-23

Given a list of SNP positions, file with annotations and directory with SNP LD files:
    - Compute the enrichment of provided set of SNPs within specified annotations
      based on annotation shifting.
    - Write the output to tab-delimited files.


Usage:
    ./goshifter.py [--snpmap FILE --proxymap FILE] --annotation FILE --permute INT --ld DIR --out FILE [--rsquared NUM --window NUM --min-shift NUM --max-shift NUM --ld-extend NUM --no-ld]

Options:
    -h, --help              Print this message and exit.
    -v, --version           Print the version and exit.

    -s, --snpmap FILE       File with SNP mappings, tab delimited, must include header: SNP, CHR, BP. Chromosomes in format chrN.
    -y, --proxymap FILE     File with proxy mappings (should include all proxies).
    -a, --annotation FILE   File with annotations, bed format. No header.
    -p, --permute INT       Number of permutations.
    -l, --ld DIR            Directory with LD files
    
    -r, --rsquared NUM      Include LD SNPs at rsquared >= NUM [default: 0.8]
    -w, --window NUM        Window size to find LD SNPs [default: 5e5]
    -n, --min-shift NUM     Minimum shift [default: False]
    -x, --max-shift NUM     Maximum shift [default: False]
    -e, --ld-extend NUM     Fixed value by which to extend LD boundries [default: False]
    -n, --no-ld             Do not include SNPs in LD [default: False]

    -o, --out FILE          Write output file.


Copyright (C) 2014  Gosia Trynka

This file is part of GoShifter package.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

BY USING THE SOFTWARE YOU ACKNOWLEDGE THAT YOU HAVE READ AND UNDERSTAND THE
TERMS OF USE BELOW. 

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

THIS SOFTWARE IS TO BE USED AS A RESEARCH TOOL ONLY. THE SOFTWARE TOOL SHALL
NOT BE USED AS A DIAGNOSTIC DECISION MAKING SYSTEM AND MUST NOT BE USED TO
MAKE A CLINICAL DIAGNOSIS OR REPLACE OR OVERRULE A LICENSED HEALTH CARE
PROFESSIONAL'S JUDGMENT OR CLINICAL DIAGNOSIS. ANY ACTION THE RESEARCHER TAKES
IN RESPONSE TO THE INFORMATION CONTAINED WITHIN IS AT THE RESEARCHER'S
DISCRETION ONLY.

"""
from __future__ import division
from datetime import datetime,date,time
from docopt import docopt
import data
import functions
import sys
import os

def warn(message):
    sys.stderr.write(message)
    sys.stderr.flush()


def invalid_arg(args, arg):
    sys.exit('ERROR: Invalid {} {}\n'.format(arg, args[arg]))


def validate_args(args):

    if not (args['--proxymap'] is None) and not os.path.exists(args['--proxymap']):
        invalid_arg(args, '--proxymap')
    elif not (args['--snpmap'] is None) and not os.path.exists(args['--snpmap']):
        invalid_arg(args, '--snpmap')


    if not os.path.exists(args['--annotation']):
        invalid_arg(args, '--annotation')

    try:
        args['--permute'] = int(args['--permute'])
    except ValueError:
        invalid_arg(args, '--permute')

    
    #if args['--ld'] != 'False':
    if not os.path.exists(args['--ld']):
        invalid_arg(args, '--ld')

    try:
        args['--rsquared'] = float(args['--rsquared'])
    except ValueError:
        invalid_arg(args, '--rsquared')

    if args['--rsquared'] <= 0 or args['--rsquared'] > 1:
        invalid_arg(args, '--rsquared')
    

    try:
        args['--window'] = float(args['--window'])
    except ValueError:
        invalid_arg(args, '--window')

    if args['--window'] <= 0:
        invalid_arg(args, '--window')
    

    if args['--min-shift'] != "False":
        try: 
            args['--min-shift'] = int(args['--min-shift'])
        except ValueError:
            invalid_arg(args, '--min-shift')
    
    if args['--max-shift'] != "False":
        try:
            args['--max-shift'] = int(args['--max-shift'])
        except:
            invalid_arg(args, '--max-shift')

    if args['--ld-extend'] != "False":
        try:
            args['--ld-extend'] = int(args['--ld-extend'])
        except:
            invalid_arg(args, '--ld-extend')

    return args


if __name__ == '__main__':
    args = docopt(__doc__,version='GoShifter 0.2')
    
    print "\n******* Analysis started", datetime.now()\
            .strftime("%A, %d. %B %Y %I:%M%p"), "*******\n"
    
    print "Running GoShifter with following parameters:\n"
    
    for arg in args:
        print "\t", arg, args[arg]
    print "\n"
    
    args = validate_args(args)


    ### read pheno mappings
    if not (args['--snpmap'] is None): 
         snpInfoChr = data.readSnpMapByChr(args['--snpmap'])
    elif not (args['--proxymap'] is None):
         snpInfoChr = data.readProxyMap(args['--proxymap'])

    #### find median annotation size to expand LD boundries
    if args['--ld-extend'] == "False":
        expand = functions.features(args['--annotation'])*2
    else:
        expand = args['--ld-extend']

    ### read peaks as interval tree
    peaksTree = functions.intervalTree(args['--annotation'])

    ### test for enrichment with peak shifting (random shift)
    functions.enrichPermRandBoundryPeakShift_tabixLd(snpInfoChr,args['--ld'],args['--rsquared'],args['--window'],expand,peaksTree,args['--min-shift'],args['--max-shift'],args['--permute'],args['--out'],args['--no-ld'])

    print "\n******* Analysis ended", datetime\
            .now().strftime("%A, %d. %B %Y %I:%M%p"), "*******\n"



