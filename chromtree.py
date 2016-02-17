"""
:File:      chromtree.py
:Author:    Kamil Slowikowski
:Updated:   October 29, 2013

This module is a simple extension of

    bx.intervals.intersection.IntervalTree

for the purpose of easily handling a separate tree for each chromosome in
a genome.

Here's a usage example:

    from chromtree import ChromTree

    tree = ChromTree()

    # A BED file with intervals on many different chromosomes.
    for line in open("file.bed"):
        chrom, begin, end, name = line.rstrip().split("\t")[:4]
        tree.insert(chrom, int(begin), int(end), name)

    # Retrieve some intervals.
    tree.find("chr1", 1000000, 1100000)

    # Output: ['54991', '254173']

"""

from bx.intervals.intersection import IntervalTree
from collections import defaultdict

class ChromTree(object):
    """A slight extension to bx.intervals.intersection.IntervalTree that can
    handle intervals on different chromosomes."""
    def __init__(self):
        self.trees = defaultdict(IntervalTree)

    def insert(self, chrom, beg, end, value):
        chrom = self._format_chrom(chrom)
        return self.trees[chrom].insert(beg, end, value)

    def find(self, chrom, beg, end):
        chrom = self._format_chrom(chrom)
        return self.trees[chrom].find(beg, end)

    def _format_chrom(self, chrom):
        """Enforce consistent chromosome names."""
        if not type(chrom) == str or not chrom.startswith('chr'):
            chrom = 'chr' + str(chrom)
        return chrom
