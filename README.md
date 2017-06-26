# GoShifter

GoShifter is a method to determine enrichment of annotations in GWAS significant loci. The manual can be found in the GoShifter_manual.pdf file. These files were taken from http://www.broadinstitute.org/mpg/goshifter/, and may have small updates and bugfixes.

Please note this code depends on bx-python: https://bitbucket.org/james_taylor/bx-python/wiki/HowToInstall
Also, this program requires tabix (http://www.htslib.org/doc/tabix.html) to be in the $PATH.

This program depends on pre-calculated LD files. Each pair of variants should be on a single line, values should be tab-separated, and the files should be indexed with Tabix. 

Format (chrA\tposA\trsIdA\tposB\trsIdB\tRsquared\tDPrime), e.g.:
chr7	143099133	rs10808026	143103481	rs56402156	0.970556	0.992515

