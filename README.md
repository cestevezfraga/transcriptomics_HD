# transcriptomics_HD
Analysis of cortical atrophy and gene expression in HD.

Code used in "Genetic topography and cortical cell loss in Huntington’s disease link development and neurodegeneration" by Estevez-Fraga et al.

The steps below can be used to investigate the transcriptomic correlates of brain maps, as far as it is in MNI space. I used the SPM T maps for differences in volume and in mean diffusivity between Huntington's patients and controls. 

I am also adding the script (a) to warp diffusion images into MNI space to perform VBM-like voxelwise analysis in SPM.

This repository contains the code to:

a) Transform MD images from native space to MNI space (dwi2md.sh) - In collaboration with Christopher S Parker, adapted from https://github.com/csparker/DWI-scripts

b) Obtain a list of genes associated with cortical charactetistics (AIBS_HD) - Developed by Andre Altmann, adapted from https://github.com/andrealtmann/AIBS_FTD 

c) Perform GO (GO) and EWCE (EWCE) in the thresholded gene lists from b) - Developed by Peter McColgan

d) Investigate the association between the expression of HTT and the expression of developmental genes using permutations from random genes (permutations) or through the spin test - Developed by Angelika Zarkali, adapted from https://github.com/AngelikaZa/YAS_HD/tree/main/Permutations

e) Create  figures summarising the information obtained with the steps above
