#!/usr/bin/env python
"""
:File: data.py
:Author: Gosia Trynka <gosia@broadinstitute.org>
:Last updated: 2014-07-23

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
import gzip
import io
import sys
import os
import re


def readSnpMapByChr(file,maf=None):
    """
    Reads in file with SNP mappings.
    Returns dict, key:chrom, value:snp, bp.
    """
    maf = maf if maf is not None else 0
    print "Filtering SNPs if MAF < {}".format(maf)
    snpMap = {}
    i = 0
    with open(file,'rU') as f:
        snpNu = 0
        header = next(f).rstrip().split("\t")
        chr_i = getCol('Chrom',header,file)
        snp_i = getCol('SNP',header,file)
        bp_i = getCol('BP',header,file)
        #check for freq in the header only if maf is specified
        if maf != 0:
            freq_i = getCol('Freq',header,file)
		
        for line in f:
            i += 1
            line = line.rstrip().split("\t")
            checkMappingsFormat(file,line,i)
            chr = line[chr_i]
            snp = line[snp_i]
            if snp == 'excluded': continue
            bp = int(line[bp_i])
            if maf != 0:
                freq = float(line[freq_i])
                if freq > 0.5:
                    freq = 1-freq
                if freq < maf: continue
            snpMap.setdefault(chr,{})
            snpMap[chr].setdefault(snp,{})
            snpMap[chr][snp]['bp'] = bp
            snpNu += 1
    print 'Read {snpNu} SNPs from {file}'.format(**locals())
    return snpMap
    
def readProxyMap(file,maf=None):
    """
    Reads in file with SNP mappings.
    Returns dict, key:chrom, value:snp, bp.
    """
    maf = maf if maf is not None else 0
    print "Filtering SNPs if MAF < {}".format(maf)
    snpMap = {}
    i = 0
    with open(file,'rU') as f:
        snpNu = 0
        header = next(f).rstrip().split("\t")
        chr_i = getCol('Chrom',header,file)
        snp_i = getCol('SNP',header,file)
        bp_i = getCol('BP',header,file)
        #check for freq in the header only if maf is specified
        if maf != 0:
            freq_i = getCol('Freq',header,file)
		
        for line in f:
            i += 1
            line = line.rstrip().split("\t")
            checkMappingsFormat(file,line,i)
            chr = line[chr_i]
            snp = line[snp_i]
            if snp == 'excluded': continue
            bp = int(line[bp_i])
            if maf != 0:
                freq = float(line[freq_i])
                if freq > 0.5:
                    freq = 1-freq
                if freq < maf: continue
            snpMap.setdefault(chr,{})
            snpMap[chr].setdefault(snp,{})
            snpMap[chr][snp]['bp'] = bp
            snpNu += 1
    print 'Read {snpNu} SNPs from {file}'.format(**locals())
    return snpMap    


def getCol(name,line,file):
    """
    Retrives the index of the column that contains queried term
    """
    try:
        col = line.index(name)
    except ValueError:
        sys.exit("File: {file} missing {name} value in the header line. STOPPING".format(**locals()))
    return col


def checkMappingsFormat(file, line, i):
    """
    Check if SNP mappings file is in the correct format.
    """
    try:
        snp = line[0]
        chr = line[1]
        bp = line[2]
        if not isSnp(snp):
            sys.exit('{file} wrong format, {snp} does not look like a SNP.'.format(**locals()))
        elif not isChr(chr):
            sys.exit("Line {i}, {file}: wrong format, {chr} does not match"\
                    " 'chrN' format. SNPs on chrom X cannot be tested."\
                    .format(**locals()))
        elif not isBp(bp):
            sys.exit("Line {i}, {file}: wrong format, {bp} does not look"\
                    " like a genomic position.".format(**locals()))
    except IndexError:
        print line
	sys.exit('{file} missing value at line {i}, is file tab delimited?'.format(**locals()))


def find(pat, text):
    foundMatch = []
    match = re.search(pat, text)
    if match:
        foundMatch = match.group()
    return foundMatch

def isHeader(snp):
	snp = find('SNP',snp)
	return snp

def isSnp(snp):
    snp = find('[\w:-]+', snp)
    return snp

def isChr(chr):
    chr = find('chr\w+',chr) ##modif to include chrX
    return chr

def isBp(bp):
    bp = find('\d+',bp)
    return bp

def isLine(line,file):
    yey = find('^chr',line)
    print yey, file

def startsChr(line):
    yay = find('^chr\w+',line)
    return yay


