# scREALTIME
A method for identifying stimulus-induced gene expression dynamics in single cells

## Dependencies
- Seurat
= Nbclust
- matrixStats
= gridExtra

## Install
Tested compatilbility with R version >=4.2.1.
Install with: 
library(devtools)\
install_github("ksheu/scResponseDynamics/scREALTIME")\
library(scREALTIME)

## Quick Start
We can use the macrophage example data provided in the 'output' folder to run scREALTIME.

### Example Use
```
macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
metadata = getMetaData(macro, stimulus = "LPS", timepoints= select_timepoints)
select_timepoints = c(0.0, 0.25, 1, 3, 8)
reconst = scREALTIME(macro = macro, metadata = metadata, timepoints = select_timepoints,
							num_archetypes = 20, data = "ISNorm",
							num_trajectories = 1000, num_sim_pts = 100,
							reduction = 'pca', stimulus = "LPS", consensus_measure = 'median',
							interpolant = 'spline', prob_method = 'distance', distance_metric = 'euclidean' ,
							varFilter = T, exp_prob = 1) 

```
scREALTIME trajectories are stored in the slot 'reconst$reconstructed_trajectories'.

### Parameters
Several parameters can be tuned by the user based on the characteristics of the dataset. 
- num_archetypes: controls the number of cell subclusters at each timepoints. A larger number will better account for the single-cell distribution at each timepoint. A smaller number will minimize the capture of outlier behavior in the resulting ensemble of trajectories. 
- num_sim_points: number of interpolated data points across the timecourse. A larger number will provide a smoother trajectory. 
- num_trajectories: number of total cells expected in the population. scREALTIME assumes that scRNAseq data is from a short timecourse, with minimal cell death or division. Too small a number may not fully explore the space allowed by the transition probability matrix that links cell archetypes across timepoints.