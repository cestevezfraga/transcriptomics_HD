# transcriptomics_HD
Analysis of cortical atrophy and gene expression in HD.

Code used in "Genetic topography and cortical cell loss in Huntington’s disease link development and neurodegeneration" by Estevez-Fraga et al.

This Repository contains the code to:

a) Transform MD images from native space to MNI space (dwi2md.sh)
b) Obtain a list of genes associated with cortical charactetistics (AIBS_HD)
c) Perform GO (GO) and EWCE (EWCE) in the thresholded gene lists from b)
d) Investigate the association between the expression of HTT and the expression of developmental genes using permutations from random genes (permutations) or through the spin test
e) Create nice figures summarising the information obtained with the steps above
