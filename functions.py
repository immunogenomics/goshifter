#!/usr/bin/env python
"""
:File: functions.py
:Author: Gosia Trynka <gosia@broadinstitute.org>
:Last updated: 2016-03-22

Copyright (C) 2014, 2015, 2016  Gosia Trynka

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
from bisect import bisect_left
from subprocess import Popen, PIPE
from chromtree import ChromTree
from bx.intervals.cluster import ClusterTree
from collections import defaultdict
import random
import copy
import sys
import os
import numpy as np
import base64
import gzip

seed=os.urandom(8)
random.seed(seed)

## to use the same seed, have to decode it back
print "\nSeed for random set to: {}\n".format(base64.b64encode(seed))

def features(fin,sep=None):
    """
    Calculates basic stats for tested annotations.
    """

    sep = sep if sep is not None else "\t" #defaults to tab delimited
    all = []
    for line in gzip.open(fin, 'r'):
        chrom, start, end = line.rstrip().split(sep)[:3]
        if chrom == "Chrom": continue
        all.append(int(end)-int(start))
    medianSize = np.median(all)
    minSize = np.min(all)
    maxSize = np.max(all)
    print "Annotation size: median = {}, min = {}, max = {}".format(medianSize,minSize,maxSize)
    return medianSize


def intervalTree(fin):
    """
    Reads gzipped bed file as interval tree
    """
    i = 0 
    tree = ChromTree()
    for line in gzip.open(fin, 'r'):
        i+= 1
        chrom, start, end = line.rstrip().split("\t")[:3]
        if chrom == "Chrom": continue
        tree.insert(chrom, int(start), int(end), (int(start),int(end)))
    print "Read {i} lines from {fin}".format(**locals())
    return tree


def enrichPermRandBoundryPeakShift_tabixLd(snpInfoChr,tabixDir,r2min,window,expand,peaksTree,minShift,maxShift,nPerm,fOut,nold):
    """
    Calculates the enrichement using random peak shifting, 
    peaks flip around if after shift they fall outside LD boundries
    """
    
    if nold:
        ldInfo = mapToLdInfo(snpInfoChr)
    else:
        ldInfo = getLdInfoTabix(snpInfoChr,tabixDir,window,r2min)
    
    snpPeakInfo, ldsnpPeakInfo = permRandBoundryPeakShift(snpInfoChr,peaksTree,ldInfo,expand,minShift,maxShift,nPerm)
    overlapPerm(snpPeakInfo,ldsnpPeakInfo,snpInfoChr,fOut,nPerm,ldInfo)


def mapToLdInfo(snpInfoChr):
    """
    If LD is not included this function creates a dictionary with index SNPs
    as its own linked SNPs.
    """
    
    ldInfo = {}
    for chrom in snpInfoChr:
        for snp in snpInfoChr[chrom]:
            bp = snpInfoChr[chrom][snp]['bp']
            ldInfo.setdefault(snp,{})
            ldInfo[snp].setdefault(snp,{})
            ldInfo[snp][snp]['bp2'] = bp
            ldInfo[snp][snp]['r2'] = 1
    return ldInfo
        

def getLdInfoTabix(snpInfoChr,tabixDir,window,r2min):
    """
    Calls tabix to retrive SNPs within window around the queried SNP. 
    Finds all the SNPs in LD at a given r2.
    """
    ldInfo = {}
    allSnps = 0

    for chrom in snpInfoChr:
        allSnps += len(snpInfoChr[chrom])
    
    i=0
    all_ld_snps = 0
    for chrom in snpInfoChr:
        for snp in snpInfoChr[chrom]:
            bp = snpInfoChr[chrom][snp]['bp']
            i += 1
            if not i % 100 or i == allSnps:
                print 'Read LD info for {i} of {allSnps} SNPs'.format(**locals())
            snpLdTabix(snp,chrom,bp,tabixDir,window,r2min,ldInfo)
            all_ld_snps += len(ldInfo[snp].keys())
    
    print "Total number of LD SNPs across {} loci is {}"\
            .format(allSnps,all_ld_snps)
    return ldInfo


def snpLdTabix(snp,chrom,bp,tabixDir,window,r2min,ldInfo):
    """
    Retrive LD info from the tabix file for a single SNP.
    """
    file = os.path.join(tabixDir,'{}.EUR.tsv.gz')
    tabixFile = file.format(chrom)
    
    ldInfo.setdefault(snp,{})
    st = bp - window
    if st < 0:
        st = 0
    en = bp + window
    
    query = "tabix {tabixFile} {chrom}:{st}-{en} | awk '$6 >= {r2min} {{print $0}}' | grep -w {snp}".format(**locals())
    proc = Popen(query,shell=True,stdout=PIPE)
    
    
    infile = 0
    for line in proc.stdout:
        fields = line.rstrip().split()
        bp1 = int(fields[1])
        snp1 = fields[2]
        bp2 = int(fields[3])
        snp2 = fields[4]
        r2 = float(fields[5])
        if snp1 == snp and r2 >= r2min:
            ldInfo[snp].setdefault(snp2,{})
            ldInfo[snp][snp2]['bp2'] = bp2
            ldInfo[snp][snp2]['r2'] = r2
            infile = 1 
        elif snp2 == snp and r2 >= r2min:
            ldInfo[snp].setdefault(snp1,{})
            ldInfo[snp][snp1]['bp2'] = bp1
            ldInfo[snp][snp1]['r2'] = r2
            infile = 1
    if not infile:
        print "{snp} was not found in {tabixFile}. "\
                "Including with no LD info.".format(**locals())
        ldInfo[snp].setdefault(snp,{})
        ldInfo[snp][snp]['bp2'] = bp
        ldInfo[snp][snp]['r2'] = 1


def countSnps(snpInfoChr):
    allSnps = 0
    for chrom in snpInfoChr:
        allSnps += len(snpInfoChr[chrom])
    return allSnps


def permRandBoundryPeakShift(snpInfoChr,peaksTree,ldInfo,expand,minShift,maxShift,nPerm):
    """
    1) Maps SNPs to shifted annotations in n permutations.
    2) Takes peaks that are mapping within the intervals by LD defined boundries. 
    3) Mentains the same number of annotations, 
       if shifted outside the boundrie the annotation flips around.
    4) Returns dict, key: iteration, SNP, value:1/0 if any of the ldSNPs 
       overlaps peak in any of the tissues in iter n.
    5) Returns dictionary with info which linked SNPs overlap annotation
    """
    
    # define LD boundries
    bounds = getLdBounds(ldInfo,snpInfoChr)

    # map peaks to LD boundries based on interval tree
    peaksToSnp = peaks2region(bounds,peaksTree,expand)
    i = 0
    snpPeakInfo = {}
    ldsnpPeakInfo = {}

    if minShift == "False" or maxShift == "False":
        print "Shifting annotations by random value within LD boundries"
    else:
        print "Shifting annotations by random value within specified "\
                "{} and {} range".format(minShift,maxShift)

    for chrom in bounds:
        for snp in bounds[chrom]:
            i += 1
            # if there are no peaks mapping to this SNP, append 0 for all iterations
            if not peaksToSnp[snp]:
                # LD SNP overlap info
                ldsnpPeakInfo.setdefault(snp,{})
                for ldsnp in ldInfo[snp]:
                    ldsnpPeakInfo[snp].setdefault(ldsnp,0)
                
                for n in range(0,nPerm+1):
                    # general overlap
                    snpPeakInfo.setdefault(n,{})            
                    snpPeakInfo[n].setdefault(snp,0)

                continue
            
            maxbp = bounds[chrom][snp]['maxbp'] + expand
            minbp = bounds[chrom][snp]['minbp'] - expand
            
            if minShift == "False" or maxShift == "False":
                #shift by random value within the size of the region
                min_shift = 0
                max_shift = maxbp-minbp
            else:
                #shift by specified value
                min_shift = minShift
                max_shift = maxShift
                
            for n in range(0,nPerm+1):
                hit = 0
                snpPeakInfo.setdefault(n,{})            
                snpPeakInfo[n].setdefault(snp,0)
                
                # for locus overlap
                ldsnpPeakInfo.setdefault(snp,{})

                if n == 0:
                    shiftVal = 0
                else:
                    shiftVal = randShift(min_shift,max_shift)

                # shift ldsnps
                for ldsnp in ldInfo[snp]:
                    ldbp = ldInfo[snp][ldsnp]['bp2'] + shiftVal
                    
                    
                    # check if shifted snp maps outside ld boundries
                    if ldbp > maxbp:
                        ldbp = ldbp - maxbp + minbp #flip around
                    elif ldbp < minbp:
                        ldbp = maxbp - (minbp - ldbp) # flip around
                    
                    # SNPs overlapping peaks in observed
                    in_peak = binarySearch(peaksToSnp[snp],ldbp)
                    if n == 0:
                        # specific SNP overlap
                        ldsnpPeakInfo[snp].setdefault(ldsnp,0)

                        if in_peak:
                            hit = 1
                            ldsnpPeakInfo[snp][ldsnp] = 1
                    # overlap in permutations
                    else:
                        if in_peak:
                            hit = 1
                            break

                snpPeakInfo[n][snp] = hit

    return snpPeakInfo, ldsnpPeakInfo


def getLdBounds(ldInfo,mappingsByChr):
    """
    For each SNP defines LD boundries based on most upstream and downstream LD SNPs
    """

    print "Defining LD boundries"
    bounds = {}
    for chrom in mappingsByChr:
        bounds.setdefault(chrom,{})
        for snp in mappingsByChr[chrom]:
            positions = []
            bounds[chrom].setdefault(snp,{})
            for ldsnp in ldInfo[snp]:
                ldbp = ldInfo[snp][ldsnp]["bp2"]
                positions.append(ldbp)
            minbp = min(positions)
            maxbp = max(positions)
            bounds[chrom][snp]['minbp'] = minbp
            bounds[chrom][snp]['maxbp'] = maxbp
    return bounds


def peaks2region(bounds,peaksTree,expand):
    """
    Assignes peaks to an interval based on bx.python interval tree.
    """

    print "Mapping annotations to LD boundries extended by {expand} bp".\
            format(**locals())
    
    peaksToSnp = {}
    for chrom in bounds:
        for snp in bounds[chrom]:
            peaksToSnp.setdefault(snp,[])
            minbp = bounds[chrom][snp]['minbp'] - expand
            maxbp = bounds[chrom][snp]['maxbp'] + expand
            peaks = peaksTree.find(chrom,minbp,maxbp)
            for peak in peaks:
                sta = peak[0]
                end = peak[1]
                # only end of peak in the bounds
                if end >=minbp and end <= maxbp and sta < minbp:
                    peaksToSnp[snp].append((minbp,end))
                # only begining of peak in the bounds
                elif sta <= maxbp and sta >= minbp and end > maxbp:
                    peaksToSnp[snp].append((sta,maxbp))
                # peak larger than bounds
                elif sta <= minbp and end >= maxbp:
                    peaksToSnp[snp].append((minbp,maxbp))
                # the entire peak is in
                elif sta >= minbp and end <= maxbp:
                    peaksToSnp[snp].append(peak)
                else:
                    sys.exit("Annotation: {sta}:{end} in region {minbp}:{maxbp}"\
                            " doesn't match any condition. Stopping!"\
                            .format(**locals()))
            peaksToSnp[snp].sort()
    return peaksToSnp


def randShift(minShift,maxShift=None):
    """
    Returns random pos/neg value between 0 and maxShift 
    """
    #min for peak shift uniform dist defaults to 0 if not specified
    minShift = minShift if minShift is not None else 0 
    x = [-1,1]
    shiftVal = random.randrange(minShift,maxShift)*random.choice(x)
    return shiftVal


def binarySearch(a, x, shift=None,lo=0,hi=None):
    """
    Binary search on the list of tuples. 
    Shift defaults to zero if not specified.
    Return the position of the tuple for which x maps between its two values.
    """
    shift = shift if shift is not None else 0 #if shift is not specified default to zero
    hi = hi if hi is not None else len(a) # hi defaults to len(a)
    x = x + shift #shift look up x if shift specified = peak shifting
    pos = bisect_left(a,(x,),lo,hi) # search list of tuples a on the 1st element
    if pos != hi and a[pos][0] == x: #position of the exact match
        return a[pos]
    elif pos-1 != hi and a[pos-1][0] < x and a[pos-1][1] >= x: # will return right index of the nearest match
        return a[pos-1]
    else:
        return False


def overlapStratPerm(snpPeakInfo,snpInfoChr,fOut,nPerm,ldInfo):
    """
    For each of n permutations, count how many SNPs overlap peaks.
    Write results to an output file. 
    """
    perm = {}
    allSnps = countSnps(snpInfoChr)

    for n in snpPeakInfo:
        count = 0 
        
        ## observed scores
        if n == 0:
            ## observed enrichment
            for snp in snpPeakInfo[n]:
                ## count total enrichment
                count += snpPeakInfo[n][snp]
                
        ## permuted scores
        else:
            for snp in snpPeakInfo[n]:
                ## count total enrichment
                count += snpPeakInfo[n][snp]
        
        perm[n] = count

    f_enrich = fOut+".nperm"+str(nPerm)+'.enrich'
    
    with open(f_enrich,'w') as f:
        f.write('nperm\tnSnpOverlap\tallSnps\tenrichment\n')
        obs = perm[0]
        betterScore = 0
        for n in perm:
            snpsOverlap = perm[n]
            enrich = round(perm[n]/allSnps,5)
            f.write('%s\t%s\t%s\t%s\n' % (n,snpsOverlap,allSnps,enrich))
            if n > 0:
                if snpsOverlap >= obs:
                    betterScore += 1

    # calculate p-value and write to log file
    pval = betterScore/nPerm
    print "p-value = {}".format(pval)
    print 'Detailed enrichment results per permutation written to {}'\
            .format(f_enrich)


def overlapPerm(snpPeakInfo,ldsnpPeakInfo,snpInfoChr,fOut,nPerm,ldInfo):
    """
    For each of n permutations, count how many SNPs overlap peaks.
    Write results to an output file. 
    """
    perm = {}
    allSnps = countSnps(snpInfoChr)
    locus_score = {}

    for n in snpPeakInfo:
        count = 0 
        
        ## observed scores
        if n == 0:
            ## observed enrichment
            for snp in snpPeakInfo[n]:
                locus_score.setdefault(snp,{})
                locus_score[snp].setdefault('overlap',0)
                locus_score[snp].setdefault('locus_score',0)
                locus_score[snp]['overlap'] = snpPeakInfo[n][snp]
            
                ## count total enrichment
                count += snpPeakInfo[n][snp]
                
        ## permuted scores
        else:
            for snp in snpPeakInfo[n]:
                locus_score[snp]['locus_score'] += snpPeakInfo[n][snp]
            
                ## count total enrichment
                count += snpPeakInfo[n][snp]
        
        perm[n] = count

    f_enrich = fOut+".nperm"+str(nPerm)+'.enrich'
    f_locus_score = fOut+".nperm"+str(nPerm)+'.locusscore'
    f_snp_score = fOut+".nperm"+str(nPerm)+'.snpoverlap'

    with open(f_enrich,'w') as f:
        f.write('nperm\tnSnpOverlap\tallSnps\tenrichment\n')
        obs = perm[0]
        betterScore = 0
        for n in perm:
            snpsOverlap = perm[n]
            enrich = round(perm[n]/allSnps,5)
            f.write('%s\t%s\t%s\t%s\n' % (n,snpsOverlap,allSnps,enrich))
            if n > 0:
                if snpsOverlap >= obs:
                    betterScore += 1

    with open (f_locus_score, 'w') as g, open(f_snp_score, 'w') as h:
        g.write('locus\toverlap\tscore\n')
        h.write('locus\tld_snp\toverlap\n')
        for snp in locus_score:
            ## if locus overlaps annotation site
            if locus_score[snp]['overlap']:
                g.write('%s\t%s\t%s\n'\
                         % (snp,locus_score[snp]['overlap'],\
                        locus_score[snp]['locus_score']/nPerm))
                for ldsnp in ldsnpPeakInfo[snp]:
                    h.write('%s\t%s\t%s\n'\
                            % (snp,ldsnp,ldsnpPeakInfo[snp][ldsnp]))
            else:
                g.write('%s\t%s\t%s\n'\
                         % (snp,'N/A','N/A'))


    # calculate p-value and write to log file
    pval = betterScore/nPerm
    print "p-value = {}".format(pval)
    print 'Detailed enrichment results per permutation written to {}'\
            .format(f_enrich)
    print 'Locus scores written to {}'\
            .format(f_locus_score)
    print 'SNP scores written to {}'\
            .format(f_snp_score)


def mergeTree(fin1):
    """
    Read bed files and merge overlapping intervals
    """
    tree = defaultdict(lambda: ClusterTree(0, 1))
    iterFileTree(fin1,tree)
    mergeIntervals = 0
    for chrom in tree:
        mergeIntervals += len(tree[chrom].getregions())
    print "There are {mergeIntervals} regions after merging overlaping"\
    " intervals".format(**locals())
    return tree


def iterFileTree(fin,trees):
    n=0 
    for line in gzip.open(fin, 'r'):
        n+=1
        chrom, beg, end = line.rstrip().split("\t")[:3]
        if chrom == "Chrom": continue
        trees[chrom].insert(int(beg), int(end), n)
    print "Read {n} lines from {fin}".format(**locals())


def merge2IntervalTree(mergeTree):
    convTree = ChromTree()
    i = 0
    for chrom, trees in mergeTree.iteritems():
        for beg, end, ns in mergeTree[chrom].getregions():
            convTree.insert(chrom,beg,end,(beg,end))
    return convTree


def enrich_shift_conditional_tabixLd(snpInfoChr,tabixDir,r2min,window,expand,a_tree,b_tree,minShift,maxShift,nPerm,fOut,nold):
    """
    Calculates the enrichment of annotation a stratifing on annotation b by random 
    peak shift, peaks flip around if shift falls outside LD boundrie 
    """
    if nold:
        ldInfo = mapToLdInfo(snpInfoChr)
    else:
        ldInfo = getLdInfoTabix(snpInfoChr,tabixDir,window,r2min)
    
    snpPeakInfo = shift_conditional(snpInfoChr,a_tree,b_tree,ldInfo,expand,minShift,maxShift,nPerm) 
    overlapStratPerm(snpPeakInfo,snpInfoChr,fOut,nPerm,ldInfo)


def shift_conditional(snpInfoChr,a_tree,b_tree,ldInfo,expand,minShift,maxShift,nPerm):
    """
    Distinguishes between enrichment being driven by one annotation (a) over
    other (b).
    Annotations a falling within b are shited within all b annotation merged
    into one segment and all a outside are shifted within single segment
    unmapped by b.
    """

    # define LD boundries
    bounds = getLdBounds(ldInfo,snpInfoChr)

    # map peaks to LD boundries
    a_peaks2Snp = peaks2region(bounds,a_tree,expand)
    b_peaks2Snp = peaks2region(bounds,b_tree,expand)

    i = 0
    # define the results dictionary
    snpPeakInfo = {}
    
    count_perm = 0
    
    if minShift == "False" or maxShift == "False":
        print "Shifting annotations by random value within LD boundries"
    else:
        print "Shifting annotations by random value within specified "\
                "{} and {} range".format(minShift,maxShift)
    
    # set results to non-overlap up front
    for chrom in bounds:
        for snp in bounds[chrom]:
            for n in range(0,nPerm+1):
                snpPeakInfo.setdefault(n,{})
                snpPeakInfo[n][snp] = 0

    for chrom in bounds:
        for snp in bounds[chrom]:
            
            # zero if no annotation a maps to this snp
            if not a_peaks2Snp[snp]:
                # test another snp
                continue
            
            # no background annotation
            elif not b_peaks2Snp[snp]:
                maxbp = bounds[chrom][snp]['maxbp'] + expand
                minbp = bounds[chrom][snp]['minbp'] - expand
            
                if minShift == "False" or maxShift == "False":
                    #shift by random value within the size of the region
                    min_shift = 1
                    max_shift = maxbp-minbp
                else:
                    #shift by specified value
                    min_shift = minShift
                    max_shift = maxShift
                
                # test overlap within LD boundries
                for n in range(0,nPerm+1):
                    # observed, no shifting
                    if n == 0:
                        shiftVal = 0
                    else:
                        # random shift defined by value between specified 
                        # minShift and maxShift
                        shiftVal = randShift(min_shift,max_shift)
                    
                    # test SNPs for overlap with a
                    for ldsnp in ldInfo[snp]:
                        ldbp = ldInfo[snp][ldsnp]['bp2'] + shiftVal
                        # check if shifted snp maps outside ld boundries
                        ldbp = check_in_bounds(ldbp,maxbp,minbp)

                        # test for peak overlap
                        in_peak = binarySearch(a_peaks2Snp[snp],ldbp)
                        if in_peak:
                            snpPeakInfo[n][snp] = 1
                            break
            else:
                maxbp = bounds[chrom][snp]['maxbp'] + expand
                minbp = bounds[chrom][snp]['minbp'] - expand
            
                # assign SNPs to B and Genome
                snps_in_b,snps_in_g = \
                        map_snps_to_annot(ldInfo[snp],b_peaks2Snp[snp])
                   
                # map a to b and genome
                a_in_b, a_in_g = \
                        find_interval_overlap(a_peaks2Snp[snp],b_peaks2Snp[snp])
                
                # make list of genomic intervals not overlapping with B
                g_peaks2Snp = \
                        non_annotation_intervals(b_peaks2Snp[snp],minbp,maxbp)
                
                # stitch and remap b intervals
                remap_b, remap_a_in_b, remap_snps_in_b, remap_b_maxbp = \
                segment_remap(b_peaks2Snp[snp],a_in_b,snps_in_b,minbp,maxbp)
                remap_b_minbp = 0
                # stitch and remap g intervals
                remap_g, remap_a_in_g, remap_snps_in_g, remap_g_maxbp = \
                segment_remap(g_peaks2Snp,a_in_g,snps_in_g,minbp,maxbp)
                remap_g_minbp = 0
                
                # do shifting
                for n in range(0,nPerm+1): 
                    # ensure there are A in B and in G for testing, and there\
                    # are SNPs in B to test
                    if remap_a_in_b and remap_a_in_g:
                        test = shift_test_overlap\
                        (n,remap_b_minbp,remap_b_maxbp,minShift,maxShift,\
                        remap_snps_in_b,remap_a_in_b)
                        
                        # next n if overlapped, else test in g
                        if test:
                            snpPeakInfo[n][snp] = test
                            continue
                    
                        # test for overlap within G
                        else:
                            test = shift_test_overlap\
                            (n,remap_g_minbp,remap_g_maxbp,minShift,maxShift,\
                            remap_snps_in_g,remap_a_in_g)
                            snpPeakInfo[n][snp] = test
                    
                    # check if there are A in B but not G
                    elif remap_a_in_b and not remap_a_in_g:
                        test = shift_test_overlap\
                        (n,remap_b_minbp,remap_b_maxbp,minShift,maxShift,\
                        remap_snps_in_b,remap_a_in_b)
                        snpPeakInfo[n][snp] = test
                    
                    # check if there are A in G but not B
                    elif remap_a_in_g and not remap_a_in_b:
                        test = shift_test_overlap\
                        (n,remap_g_minbp,remap_g_maxbp,minShift,maxShift,\
                        remap_snps_in_g,remap_a_in_g)
                        snpPeakInfo[n][snp] = test

                    else:
                        sys.exit("I missed something. No A annotations"\
                                " for SNP {}".fomrat(snp))

    return snpPeakInfo


def map_snps_to_annot(ldInfo,annot_list):
    """
    Maps SNPs to annotation if they overlap it, otherwise to Genome (fragments
    unmapped by annotation B).
    Returns two lists with bp of all tested SNPs.
    """
    snps_in_b = []
    snps_in_g = []
    for ldsnp in ldInfo:
        ldbp = ldInfo[ldsnp]['bp2']
        in_peak = binarySearch(annot_list,ldbp)
        if in_peak:
            snps_in_b.append(ldbp)
        else:
            snps_in_g.append(ldbp)

    snps_in_b.sort()
    snps_in_g.sort()
    return snps_in_b, snps_in_g


def find_interval_overlap(a_annot_list,b_annot_list):
    """
    Finds annotations A overlapping or mapping outside B.
    Returns two lists of tuples.
    """
    a_annot_list_cp = copy.deepcopy(a_annot_list)
    a_in_b = []
    a_in_g = []
    for b in b_annot_list:
        b_st = b[0]
        b_en = b[1]
        # this will be modified through iterations
        a_annot_list_mod = copy.deepcopy(a_annot_list_cp)
        for a in a_annot_list_mod:
            a_st = a[0]
            a_en = a[1]
            # maps outside interval b
            if a_st >= b_en:
                continue
            
            # entirely maps within interval
            elif a_st >= b_st and a_en <= b_en:
                a_in_b.append((a_st,a_en))
                a_annot_list_cp.remove(a)
            
            # end of a overlaps interval b
            elif a_st < b_st and a_en <= b_en and a_en >= b_st:
                a_in_b.append((b_st,a_en))
                # append the overhang
                a_annot_list_cp.remove(a)
                a_annot_list_cp.append((a_st,b_st))
            
            # begining of a overlaps interval b
            elif a_st >= b_st and a_st <= b_en and a_en > b_en:
                a_in_b.append((a_st,b_en))
                # append the overhang
                a_annot_list_cp.remove(a)
                a_annot_list_cp.append((b_en,a_en))
            
            # a extends interval b on both ends
            elif a_st < b_st and a_en > b_en:
                a_in_b.append((b_st,b_en))
                # append the overhang
                a_annot_list_cp.remove(a)
                a_annot_list_cp.append((a_st,b_st))
                a_annot_list_cp.append((b_en,a_en))
    
    for a in a_annot_list_cp:
        a_in_g.append(a)
    a_in_b.sort()
    a_in_g.sort()
    return a_in_b, a_in_g


def non_annotation_intervals(annot_list,minbp,maxbp):
    """
    Finds genomic regions without annotations and returns them as a list of
    tuples.
    """
    interval_list = []
    all_annot = len(annot_list)
    # only one annotation
    if all_annot == 1:
        annot = annot_list[0]
        annot_st = annot[0]
        annot_en = annot[1]
        # annot spans the whole region
        if annot_st == minbp and annot_en == maxbp:
            interval_list = []
        # end of annot in the region
        elif annot_st == minbp and annot_en < maxbp:
            interval_list.append((annot_en,maxbp))
        # start of annot in the region
        elif annot_st > minbp and annot_en == maxbp:
            interval_list.append((minbp,annot_en))
        # whole annot within the region
        elif annot_st > minbp and annot_en < maxbp:
            interval_list.append((minbp,annot_st))
            interval_list.append((annot_en,maxbp))
        # catch other situations
        else:
            sys.exit('Region limits: {minbp},{maxbp}. Annotation limits:'\
                    ' {annot_st},{annot_en}. Something weird! STOPPING\n'\
                    .format(**locals()))
    # more than one annotation
    else:
        i = 0
        for annot in annot_list:
            i += 1
            annot_st = annot[0]
            annot_en = annot[1]
            # if it's the last annotation
            if i == 1:
                # annotation starts at the region boundry
                if annot_st == minbp:
                    st = annot_en
                    continue
                else:
                    st = annot_en
                    interval_list.append((minbp,annot_st))
            elif i == all_annot:
                # if end is the same as end of the region
                if annot_en == maxbp:
                    interval_list.append((st,annot_st))
                    break
                else:
                    interval_list.append((st,annot_st))
                    interval_list.append((annot_en,maxbp))
            # annot within region
            else:
                interval_list.append((st,annot_st))
                st = annot_en
    interval_list.sort()
    return interval_list


def intervals_sum_length(intervals_list):
    length = 0
    for a in intervals_list:
        sta = a[0]
        end = a[1]
        length += end - sta
    return length

def segment_remap(intervals_list,a_in_intervals,snps_in_intervals,minbp,maxbp):
    """
    Stiches all intervals into one region.  
    SNPs and annotation A within intervals are remapped to new
    coordinates.
    """
    remap_intervals_list = []
    remap_a_in_intervals = []
    remap_snps_in_intervals = []
    remap_maxbp = 0

    all_intervals = len(intervals_list)
    all_a = len(a_in_intervals)
    all_snps = len(snps_in_intervals)
    i = 0
    
    for i in range(0,all_intervals):

        if i == 0:
            sta = intervals_list[0][0]
            end = intervals_list[0][1]
            delta = sta
            remap_sta = sta - delta
            remap_end = end - delta
            remap_intervals_list.append((remap_sta,remap_end))
            
            # remap a annotations falling into this interval
            for a in a_in_intervals:
                a_sta = a[0]
                a_end = a[1]
                # don't iterate if start maps beyon the interval
                if a_sta > end:
                    break
                elif a_sta >= sta and a_end <= end:
                    remap_a_sta = a_sta - delta
                    remap_a_end = a_end - delta
                    remap_a_in_intervals.append((remap_a_sta,remap_a_end))
            
            # remap SNPs falling into this interval
            for bp in snps_in_intervals:
                if bp > end:
                    break
                elif bp >= sta and bp <= end:
                    remap_bp = bp - delta
                    remap_snps_in_intervals.append(remap_bp)
            if all_intervals == 1:
                remap_maxbp = remap_end

        else:
            sta = intervals_list[i][0]
            delta = delta + sta - end
            end = intervals_list[i][1]
            remap_sta = sta-delta
            remap_end = end-delta
            remap_intervals_list.append((remap_sta,remap_end))

            # remap a annotations falling into this interval
            for a in a_in_intervals:
                a_sta = a[0]
                a_end = a[1]
                # don't iterate if start maps beyond the interval
                if a_sta > end:
                    break
                elif a_sta >= sta and a_end <= end:
                    remap_a_sta = a_sta - delta
                    remap_a_end = a_end - delta
                    remap_a_in_intervals.append((remap_a_sta,remap_a_end))
            
            # remap SNPs falling into this interval
            for bp in snps_in_intervals:
                if bp > end:
                    break
                elif bp >= sta and bp <= end:
                    remap_bp = bp - delta
                    remap_snps_in_intervals.append(remap_bp)
        
            if i == all_intervals-1:
                remap_maxbp = remap_end
    
    sum_a = intervals_sum_length(a_in_intervals)
    remap_sum_a = intervals_sum_length(remap_a_in_intervals)
    sum_intervals = intervals_sum_length(intervals_list)
    remap_sum_intervals = intervals_sum_length(remap_intervals_list)
    
    len_remap_intervals = len(remap_intervals_list)
    assert all_intervals == len(remap_intervals_list),\
            'Number of intervals after remaping {len_remap_intervals}'\
            ' is different than in the input {all_intervals}!'\
            ' Stopping!!'.format(**locals())
    len_remap_a = len(remap_a_in_intervals)
    if all_a != len(remap_a_in_intervals):
        print minbp,maxbp
        print "intervals:",intervals_list
        print "remapped intervals:", remap_intervals_list
        print "a in intervals:",a_in_intervals
        print "a in remapped intervals:", remap_a_in_intervals
    assert all_a == len(remap_a_in_intervals),\
            'Number of A annotations after remaping {len_remap_a}'\
            ' is different than in the input {all_a}!'\
            ' Stopping!!'.format(**locals())
    len_remap_snps = len(remap_snps_in_intervals)
    assert all_snps == len(remap_snps_in_intervals),\
            'Number of SNPs after remaping {len_remap_snps} is different than'\
            ' in the input {all_snps}! Stopping!!'.format(**locals())
    assert remap_sum_a == sum_a,\
            'Length of intervals A after remaping {remap_sum_a} is different'\
            ' than in the input {sum_a}! Stopping!! A: {a_in_intervals}'\
            ' A_remap: {remap_a_in_intervals}'.format(**locals())
    assert remap_sum_intervals == sum_intervals,\
            'Length of intervals A after remaping {remap_sum_intervals}'\
            ' is different than in the input {sum_intervals}.'\
            ' Stopping!!'.format(**locals())

    return remap_intervals_list,remap_a_in_intervals,\
        remap_snps_in_intervals,remap_maxbp


def shift_test_overlap(n,minbp,maxbp,minShift,maxShift,snps,a_list):
    
    test = 0
    # if there are no SNPs to test
    if not snps:
        return test
    else:
        # define shifting value
        if n == 0:
            # observed
            shiftVal = 0
        else:
            # no limits on shifting values
            if minShift == "False" or maxShift == "False":
                min_shift = minbp
                max_shift = maxbp
            else:
                #shift within specified value
                min_shift = minShift
                max_shift = maxShift
            shiftVal = randShift(min_shift,max_shift)
    
    for ldbp in snps:
        ldbp_shift = ldbp + shiftVal
        ldbp_shift = \
        check_in_bounds(ldbp_shift,maxbp,minbp)
        in_peak = binarySearch(a_list,ldbp_shift)
        if in_peak:
            test = 1
            break

    return test


def check_in_bounds(ldbp,maxbp,minbp):
    if ldbp > maxbp:
        ldbp = ldbp - maxbp + minbp #flip around
    elif ldbp < minbp:
        ldbp = maxbp - (minbp - ldbp) # flip around
    return ldbp


