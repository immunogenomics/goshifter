# GoShifter

GoShifter is a method to determine enrichment of annotations in GWAS significant loci. The manual can be found in the GoShifter_manual.pdf file. These files were taken from http://www.broadinstitute.org/mpg/goshifter/, and may have small updates and bugfixes.

Please note this code depends on bx-python: https://bitbucket.org/james_taylor/bx-python/wiki/HowToInstall
Also, this program requires tabix (http://www.htslib.org/doc/tabix.html) to be in the $PATH.
Furthermore, this program requires LD files that can be obtained here: https://data.broadinstitute.org/srlab/BEAGLE/1kG-beagle-release3/pairwise_ld/r2_ge_0.8/

If for any reason the LD files are not available anymore, use the following structure (tab-separated):
chromosome1 pos1 rs1 pos2 rs2 r-squared dprime

for example:
chr7	143099133	rs10808026	143103481	rs56402156	0.970556	0.992515

