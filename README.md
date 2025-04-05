# scResponseDynamics
Stimulus-induced gene expression dynamics in single macrophage cells (Sheu et al, 2024)

<img src="https://github.com/KSheu/scResponseDynamics/blob/main/GA_MolCell2024.tiff" width="350" height="350">

## In Brief
Because RNA measurements are cell destructive, it is unclear how variable stimulus-induced dynamic gene expression trajectories (GETs) are. Sheu et al. developed a method to assess single-cell gene expression trajectories (scGETs) in macrophages responding to stimuli and found scGETs to be much more informative than any single timepoint measurement.

## Outline
See installation instructions for the scREAL-TIME R package under the scREALTIME subfolder. \
The run_compile.R file contains code for the analysis in Figures 1-7. \
Figure 1: simulation of scGETs using mathematical models \
Figure 2: evaluation of imputed scGETs against simulated ground truth \
Figure 3: plotting of scRNAseq data from stimulated macrophages and imputation of scGETs\
Figure 4: calculating trajectory feature values and mutual information analysis \
Figure 5: analysis of gene-gene correlations across trajectory features for single cells \
Figure 6: comparison of single-cell response dynamics for different macrophage polarization states \
Figure 7: identification of macrophage polarization state based on single-cell response dynamics 



## Citation
Sheu, Katherine M., Aditya Pimplaskar, and Alexander Hoffmann. “Single-Cell Stimulus-Response Gene Expression Trajectories Reveal the Stimulus Specificities of Dynamic Responses by Single Macrophages.” Molecular Cell 84, no. 21 (November 2024): 4095-4110.e6. https://doi.org/10.1016/j.molcel.2024.09.023.
