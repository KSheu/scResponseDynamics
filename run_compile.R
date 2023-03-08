#scResponseDynamics_script: imputation of single-cell response trajectories

source('trajectory4.R')

library(data.table);library(reshape2);library(pheatmap);
library(RColorBrewer);library(ggpubr);library(ggplot2);
library(ksheu.library1);library(SLEMI)library(Seurat);library(ggpubr);
setwd("F://scRNAseq_macro/scRNAseq_macro/")

###############################################################################
# Figure 1: plot data ----
#plot violins, group by timepoint----
if(1){
  macro = readRDS("output/macrophage_M0_rep2only_500genes_DBEC.rds")
  gene = "Tnfsf8"      # "Tnfaip3","Nfkbiz", "Ptges"
  ymax = 10
  pt.size = 0.01
  p1=VlnPlot(object = subset(macro, subset= stimulus=="LPS"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#00BA38",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("LPS")
  p2=VlnPlot(object = subset(macro, subset= stimulus=="PIC"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#619CFF",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("PIC")
  p3=VlnPlot(object = subset(macro, subset= stimulus=="IFNb"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#B79F00",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("IFNb")
  p4=VlnPlot(object = subset(macro, subset= stimulus=="P3CSK"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#00BFC4",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("P3C")
  p5=VlnPlot(object = subset(macro, subset= stimulus=="CpG"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#F8766D",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("CpG")
  p6=VlnPlot(object = subset(macro, subset= stimulus=="TNF"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#F564E3",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("TNF")
  p.g1 = (p1|p2|p3)/(p4|p5|p6)
  p.g1
  
  macro = readRDS("output/macrophage_M1_IFNg_500genes_DBEC.rds")
  macro = subset(macro, subset = timept!= "x24hr")
  p1=VlnPlot(object = subset(macro, subset= stimulus=="LPS"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#00BA38",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("LPS")
  p2=VlnPlot(object = subset(macro, subset= stimulus=="PIC"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#619CFF",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("PIC")
  p3=VlnPlot(object = subset(macro, subset= stimulus=="IFNb"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#B79F00",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("IFNb")
  p4=VlnPlot(object = subset(macro, subset= stimulus=="P3CSK"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#00BFC4",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("P3C")
  p5=VlnPlot(object = subset(macro, subset= stimulus=="CpG"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#F8766D",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("CpG")
  p6=VlnPlot(object = subset(macro, subset= stimulus=="TNF"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#F564E3",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("TNF")
  
  p.g2 = (p1|p2|p3)/(p4|p5|p6) 
  p.g2
  
  
  macro = readRDS("output/macrophage_M2_IL4_gt80_500genes_DBEC.rds")
  p1=VlnPlot(object = subset(macro, subset= stimulus=="LPS"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#00BA38",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("LPS")
  p2=VlnPlot(object = subset(macro, subset= stimulus=="PIC"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#619CFF",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("PIC")
  p3=VlnPlot(object = subset(macro, subset= stimulus=="IFNb"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#B79F00",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("IFNb")
  p4=VlnPlot(object = subset(macro, subset= stimulus=="P3CSK"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#00BFC4",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("P3C")
  p5=VlnPlot(object = subset(macro, subset= stimulus=="CpG"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#F8766D",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("CpG")
  p6=VlnPlot(object = subset(macro, subset= stimulus=="TNF"|timept=="0.0hr"), pt.size = pt.size, cols = rep("#F564E3",7),
             features = c(gene), group.by = "timept", assay = "ISnorm" , y.max = ymax) +
    stat_summary(fun.y = median, geom='point', size = 1, colour = "green") +
    theme(legend.position = "None")+ylab(NULL)+xlab(NULL)+ggtitle("TNF")
  
  p.g3 = (p1|p2|p3)/(p4|p5|p6)
  p.g3
  p.g1|p.g2|p.g3
  
  p.g1|p.g2|p.g3 +ylab(gene)+ theme(plot.margin = unit(c(0,300,0,0), "pt"))
}



###############################################################################
# Figure 2: impute all trajectories - March 2022 redo----
# Note2self: from F://BACKUP.../Projects.../trajectory_method/tensor_trajectory.R

# for M0 macrophages-----
reduction_name = 'pca'
select_timepoints = c(0.0, 0.25,1,3,8)
macro = readRDS("./output/macrophage_M0_rep2only_500genes_DBEC.rds")
macro = FindVariableFeatures(macro,assay = "ISnorm")
macro@active.assay = "ISnorm"
macro = ScaleData(macro)
macro[["ISnorm"]]@scale.data = macro[["ISnorm"]]@data
macro <- RunPCA(macro, assay = "ISnorm", scale =F, center = F)
macro <- RunNMF(macro, assay = "ISnorm", scale =F, center = F)
macro <- RunICA(macro, assay = "ISnorm")
for (n in c("LPS", "TNF","P3CSK","CpG", "PIC","IFNb")){
  metadata = getMetaData(macro, stimulus = n, timepoints= select_timepoints)
  # PCAPlot(macro, dims = c(1,2), cells = as.vector(rownames(metadata)), group.by = 'timept', shuffle = TRUE, label = TRUE, label.box = TRUE) + ggtitle('PCA with PC 2 and 3')
  num_archetypes = 20
  
  interpl = "spline"
  # undebug(getTrajectory4)
  reconst = getTrajectory4(macro = macro, metadata = metadata, num_archetypes = num_archetypes, data = "ISNorm",
                           timepoints = select_timepoints, num_trajectories = 1000, num_sim_pts = 100,
                           reduction = reduction_name, stimulus = n, consensus_measure = 'median',
                           interpolant = interpl, prob_method = 'distance', distance_metric = 'euclidean' ,
                           varFilter = T, exp_prob = 1) # saying empty clusters -- insufficient data??
  
  # saveRDS(reconst, paste0("./trajectory/reconstructed_M0_rep2only_500genes_",reduction_name,"onallstims_",n, "_",interpl,"_k_", num_archetypes,".rds"))
}
# reconst = readRDS(paste0("./trajectory/reconstructed_M0_rep2only_500genes_",reduction_name,"_",n, "_",interpl,"_k_", num_archetypes,".rds"))

# for M1(IFNg) macrophages-----
select_timepoints = c(0.0, 0.5, 1,3,5,8)
macro = readRDS("./output/macrophage_M1_IFNg_500genes_DBEC.rds")
macro = subset(macro, timept != "x24hr")
macro = FindVariableFeatures(macro,assay = "ISnorm")
macro@active.assay = "ISnorm"
macro = ScaleData(macro)
macro[["ISnorm"]]@scale.data = macro[["ISnorm"]]@data
macro <- RunPCA(macro, assay = "ISnorm")
for (n in c("LPS", "TNF","P3CSK","CpG", "PIC","IFNb")){
  metadata = getMetaData(macro, stimulus = n, timepoints= select_timepoints)
  # PCAPlot(macro, dims = c(2,3), cells = as.vector(rownames(metadata)), group.by = 'timept', shuffle = TRUE, label = TRUE, label.box = TRUE) + ggtitle('PCA with PC 2 and 3')
  num_archetypes = 20
  
  interpl = "spline"
  # undebug(getTrajectory4)
  reconst = getTrajectory4(macro = macro, metadata = metadata, num_archetypes = num_archetypes, data = "ISNorm",
                           timepoints = select_timepoints, num_trajectories = 1000, num_sim_pts = 100,
                           reduction = reduction_name, stimulus = n, consensus_measure = 'median',
                           interpolant = interpl, prob_method = 'distance', distance_metric = 'euclidean' ,
                           varFilter = T, exp_prob = 1) # saying empty clusters -- insufficient data??
  
  saveRDS(reconst, paste0("./trajectory/reconstructed_M1_IFNg_500genes_",reduction_name,"onallstims_",n, "_",interpl,"_k_", num_archetypes,".rds"))
}

#for M2(IL4) macrophages------
select_timepoints = c(0.0, 0.5,1,3,5,8)
macro = readRDS("./output/macrophage_M2_IL4_gt80_500genes_DBEC.rds")
macro = FindVariableFeatures(macro,assay = "ISnorm")
macro@active.assay = "ISnorm"
macro = ScaleData(macro)
macro[["ISnorm"]]@scale.data = macro[["ISnorm"]]@data
macro <- RunPCA(macro, assay = "ISnorm", scale =F, center = F)
for (n in c("LPS", "TNF","P3CSK","CpG", "PIC","IFNb")){
  metadata = getMetaData(macro, stimulus = n, timepoints= select_timepoints)
  # PCAPlot(macro, dims = c(2,3), cells = as.vector(rownames(metadata)), group.by = 'timept', shuffle = TRUE, label = TRUE, label.box = TRUE) + ggtitle('PCA with PC 2 and 3')
  num_archetypes = 20
  
  interpl = "spline"
  # undebug(getTrajectory4)
  reconst = getTrajectory4(macro = macro, metadata = metadata, num_archetypes = num_archetypes, data = "ISNorm",
                           timepoints = select_timepoints, num_trajectories = 1000, num_sim_pts = 100,
                           reduction = reduction_name, stimulus = n, consensus_measure = 'median',
                           interpolant = interpl, prob_method = 'distance', distance_metric = 'euclidean' ,
                           varFilter = T, exp_prob = 1) # saying empty clusters -- insufficient data??
  
  saveRDS(reconst, paste0("./trajectory/reconstructed_M2_IL4_gt80_500genes_",reduction_name,"onallstims_",n, "_",interpl,"_k_", num_archetypes,".rds"))
}

####################################################### plot the trajectories-----
#plot scaled across all stimuli, comment/uncomment for M0, M1, M2----
reduction_name = "pca"
interpl = "spline"
num_archetypes = 20
mat_allstims <- list()
num_sim_pts = 100+3
num_paths_sum = 0
labels = c()
labels2 = c()
gc();gc();gc()
for (stim in c("LPS", "TNF","P3CSK","CpG", "PIC","IFNb")){
  print(stim)
  # reconstructed_pc = readRDS(paste0("./trajectory/reconstructed_",reduction_name,"_",stim,"_k_", num_archetypes))
  reconstructed_pc = readRDS(paste0("./trajectory/reconstructed_M0_rep2only_500genes_",reduction_name,"onallstims_",stim, "_",interpl,"_k_", num_archetypes,".rds"))
  # reconstructed_pc = readRDS(paste0("./trajectory/reconstructed_M1_IFNg_500genes_",reduction_name,"onallstims_",stim, "_",interpl,"_k_", num_archetypes,".rds"))
  # reconstructed_pc = readRDS(paste0("./trajectory/reconstructed_M2_IL4_gt80_500genes_",reduction_name,"onallstims_",stim, "_",interpl,"_k_", num_archetypes,".rds"))
  
  num_paths = length(unique(reconstructed_pc$reconstructed_trajectories$path))
  num_timepts = length(unique(reconstructed_pc$reconstructed_trajectories$time))
  
  reconstructed_pc.traj = (reconstructed_pc$reconstructed_trajectories)
  reconstructed_pc.traj$stimulus = stim
  mat_allstims <- rbind(mat_allstims, reconstructed_pc.traj)
  
  num_paths_sum = num_paths_sum + num_paths
  labels = c(labels, rep(stim, num_paths))
}

mat.numbers = mat_allstims[,!grepl("time|path|stimulus", colnames(mat_allstims))]
mat.meta = mat_allstims[,grepl("time|path|stimulus", colnames(mat_allstims))]
if (1){ #rescale 0-1 or not----
  mat.numbers = apply(mat.numbers, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X))) #rescale each gene column 0-1 over all stims
  # mat.numbers = cbind(mat.numbers, mat.meta)
} 
if (0){ #scale to max and set negative to 0----
  mat.numbers = apply(mat.numbers, MARGIN = 2, FUN = function(X) (X /max(X))) #rescale each gene column to max over all stims
  mat.numbers[mat.numbers <0] = 0
} 

#######################################################################################
#plot scaled across all stimuli and mac types for comparing polarization states-------
if(0){ #either scale each mac type individually above, or together  below
if(1){
  stim="LPS"
  reconstructed_pc = readRDS(paste0("./trajectory/reconstructed_","M0_rep2only","_500genes_",reduction_name,"onallstims_",stim, "_",interpl,"_k_", num_archetypes,".rds"))
  genes1 = (colnames(reconstructed_pc$reconstructed_trajectories))
  time1 = unique(reconstructed_pc$reconstructed_trajectories$time)
  reconstructed_pc = readRDS(paste0("./trajectory/reconstructed_","M1_IFNg","_500genes_",reduction_name,"onallstims_",stim, "_",interpl,"_k_", num_archetypes,".rds"))
  genes2 = (colnames(reconstructed_pc$reconstructed_trajectories))
  time2 = unique(reconstructed_pc$reconstructed_trajectories$time)
  reconstructed_pc = readRDS(paste0("./trajectory/reconstructed_","M2_IL4_gt80","_500genes_",reduction_name,"onallstims_",stim, "_",interpl,"_k_", num_archetypes,".rds"))
  genes3 = (colnames(reconstructed_pc$reconstructed_trajectories))
  time3 = unique(reconstructed_pc$reconstructed_trajectories$time)
}
geneset_M0M1M2 = intersect_all(genes1,genes2,genes3)
timeset_M0M1M2 = intersect_all(time1,time2,time3)

gc();gc();gc()
for (type in c("M0_rep2only","M1_IFNg","M2_IL4_gt80")){
  for (stim in c("LPS", "TNF","P3CSK","CpG", "PIC","IFNb")){
    print(type)
    print(stim)
    
    reconstructed_pc = readRDS(paste0("./trajectory/reconstructed_",type,"_500genes_",reduction_name,"onallstims_",stim, "_",interpl,"_k_", num_archetypes,".rds"))
    
    num_paths = length(unique(reconstructed_pc$reconstructed_trajectories$path))
    num_timepts = length(unique(timeset_M0M1M2))
    
    reconstructed_pc.traj = (reconstructed_pc$reconstructed_trajectories)
    reconstructed_pc.traj$stimulus = stim
    reconstructed_pc.traj$type = type
    mat_allstims <- rbind(mat_allstims, reconstructed_pc.traj[reconstructed_pc.traj$time %in% c(timeset_M0M1M2), 
                                                              c(geneset_M0M1M2,"stimulus","type")])
    
    num_paths_sum = num_paths_sum + num_paths
    labels = c(labels, rep(stim, num_paths))
    labels2 = c(labels2, rep(type, num_paths))
  }
}
mat.numbers = mat_allstims[,!grepl("time|path|stimulus|type", colnames(mat_allstims))]
mat.meta = mat_allstims[,grepl("time|path|stimulus|type", colnames(mat_allstims))]
if (1){ #rescale 0-1 or not----
  mat.numbers = apply(mat.numbers, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X))) #rescale each gene column 0-1 over all stims
  # mat.numbers = cbind(mat.numbers, mat.meta)
} 
gc();gc();gc()
}

# convert to tensor object----
dims <- c( num_timepts, num_paths_sum, ncol(reconstructed_pc$reconstructed_trajectories)-2)
dims <- c( num_timepts, num_paths_sum, length(geneset_M0M1M2)-2 )
arr <- array(as.numeric(unlist(mat.numbers)), dim = dims)
library(rTensor)
A = as.tensor(arr)
A@modes
# saveRDS(A, "~/ksheu/tensor_M0_rep2only_allstims_103x5996x497_k20.rds")
# saveRDS(A, "~/ksheu/tensor_M1_IFNg_allstims_104x6000x495_k20.rds")
# saveRDS(A, "~/ksheu/tensor_M2_IL4_allstims_104x6000x498_k20.rds")

# A = readRDS("~/ksheu/tensor_M0_rep2only_allstims_103x5996x497_k20.rds")
# A = readRDS("~/ksheu/tensor_M1_IFNg_allstims_104x6000x495_k20.rds")
# A = readRDS("~/ksheu/tensor_M2_IL4_allstims_104x6000x498_k20.rds")

####################################
## TCA Tucker decomposition-----
# A = readRDS("~/ksheu/tensor_M0_rep2only_allstims_103x5996x497_k20.rds")
library(rTensor)
A@modes;A@num_modes;dim(A)
num_timepts = A@modes[1]
num_archetypes = 20
reduction_name="pca"
# tucker_decomp <- rTensor::tucker(A, ranks = c(num_timepts, 10, 10))
str(tucker_decomp)
timepoints = c(0,0.25,1,3,8)
# saveRDS(file = paste0("~/ksheu/tensor.tuckerdecomp_reconstructed_",reduction_name, "_k_" , num_archetypes, ".rds"), tucker_decomp)
tucker_decomp = readRDS(paste0("~/ksheu/tensor.tuckerdecomp_reconstructed_",reduction_name, "_k_" , num_archetypes, ".rds"))

#plotting the TCA output------
mat.meta.agg= mat.meta[!duplicated(mat.meta[,-1]),]
U.1 = data.frame(tucker_decomp$U[[1]]);
rownames(U.1) = 1:num_timepts
p1=ggplot(U.1, aes(as.numeric(time1), X1*-1))+geom_point()+geom_path()+theme_bw(base_size = 14)
p2=ggplot(U.1, aes(as.numeric(time1), X2))+geom_point()+geom_path()+theme_bw(base_size = 14)
p3=ggplot(U.1, aes(as.numeric(time1), X3))+geom_point()+geom_path()+theme_bw(base_size = 14)
p4=ggplot(U.1, aes(as.numeric(time1), X4*-1))+geom_point()+geom_path()+theme_bw(base_size = 14)
p5=ggplot(U.1, aes(as.numeric(time1), X5*-1))+geom_point()+geom_path()+theme_bw(base_size = 14)
p1|p2|p3|p4|p5

p6=ggplot(U.1, aes(as.numeric(time1), X6*-1))+geom_point()+geom_path()+theme_bw(base_size = 14)
p7=ggplot(U.1, aes(as.numeric(time1), X7))+geom_point()+geom_path()+theme_bw(base_size = 14)
p8=ggplot(U.1, aes(as.numeric(time1), X8))+geom_point()+geom_path()+theme_bw(base_size = 14)
p9=ggplot(U.1, aes(as.numeric(time1), X9*-1))+geom_point()+geom_path()+theme_bw(base_size = 14)
p10=ggplot(U.1, aes(as.numeric(time1), X10*-1))+geom_point()+geom_path()+theme_bw(base_size = 14)
p6|p7|p8|p9|p10

U.2 = data.frame(tucker_decomp$U[[2]]);
U.2 = cbind(stimulus=mat.meta.agg$stimulus, U.2)

U.3 = data.frame(tucker_decomp$U[[3]]);
rownames(U.3) = make.unique(colnames(reconstructed_pc$reconstructed_trajectories)[1:A@modes[3]])

library(rgl);library(plot3Drgl);
colors_list = c(CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3")
with(U.2, scatter3D(x=X2, y=X3, z=X4*-1, ticktype = "detailed",nticks=4, 
                    colvar = as.integer(as.factor(U.2$stimulus)), axis.scales=T, axis.ticks = T,
                    col = colors_list, bty = "b2", pch = 19,  alpha = 0.1, cex = .75, #type = "h", 
                    xlim = c(-0.03,0.02), ylim =  c(-0.03,0.02), zlim = c(-0.03,0.03),
                    surface=T,ellipsoid=T, level=0.8, ellipsoid.alpha=0.1,groups = stimulus, surface.col = colors_list,
                    xlab = "X2", ylab="X3",zlab="X4",
                    phi = 25, theta = 45))
with(U.2, scatter3D(x=X1*-1, y=X2, z=X3, ticktype = "detailed", nticks=4,
                    colvar = as.integer(as.factor(U.2$stimulus)), axis.scales=T, axis.ticks = T,
                    col = colors_list, bty = "b2", pch = 19,  alpha = 0.1, cex = .75,# type = "h", 
                    # xlim = c(-0.03,0.02), ylim =  c(-0.03,0.02), zlim = c(-0.03,0.03),
                    surface=T,ellipsoid=T, level=0.8, ellipsoid.alpha=0.1,groups = stimulus, surface.col = colors_list,
                    xlab = "X1", ylab="X2",zlab="X3",
                    phi = 25, theta = 45))
plotrgl()

###############################################################################
#plot trajectories----
mat.numbers2 = cbind(mat.numbers, mat.meta)
mat.numbers2$path_stim = paste0(mat.numbers2$path,"_", mat.numbers2$stimulus)
mat.numbers2$path_stim = paste0(mat.numbers2$path,"_", mat.numbers2$stimulus,"_",mat.numbers2$type)
ggplot(mat.numbers2[grepl("LPS",mat.numbers2$stimulus)&grepl("",mat.numbers2$type),], aes(time,Nfkbia, group = path_stim)) + ylim(0,1)+
  # geom_vline(xintercept = c(0,0.25,1,3,8), linetype="dotted")+ #theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  geom_vline(xintercept = c(0,0.5,1,3,5,8), linetype="dotted")+
  geom_line(aes(group = as.factor(path_stim), color = type), alpha = 0.05)+
  theme_bw(base_size = 14)+theme(legend.position = "none")
#single gene------
gene = "Tnf"
# mat.numbers2=trajectories_M2
mat.numbers2.dcast = dcast(mat.numbers2, stimulus+path~time, value.var = gene)
rownames(mat.numbers2.dcast) = paste0(mat.numbers2.dcast$path,"_", mat.numbers2.dcast$stimulus)
dynamics = dynamics_M0[grepl(paste0(gene,"$"), dynamics_M0$gene),]
dynamics = dynamics_M1[grepl(paste0(gene,"$"), dynamics_M1$gene),]
dynamics = dynamics_M2[grepl(paste0(gene,"$"), dynamics_M2$gene),]
dynamics = dynamics[order(dynamics$stimulus),]
mat.numbers2.dcast = cbind(dynamics ,mat.numbers2.dcast)
colors_list = list(stimulus = c(CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3"))
annot.frame = mat.numbers2.dcast[,c(4,5,6,8,13),drop=F] #c(13,4,5,6,8,9,10)
mat.numbers2.dcast2 = mat.numbers2.dcast[grepl("", mat.numbers2.dcast$stimulus),]

# pheatmap(as.matrix(mat.numbers2.dcast2[order(mat.numbers2.dcast2$peak_amp),-c(1:14)]), 
pheatmap(as.matrix(mat.numbers2.dcast2[order(mat.numbers2.dcast2$peak_amp, decreasing = T),-c(1:14)]),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         breaks = c(0,seq(0.01,0.99,length=100),1),
         cluster_cols = F, cluster_rows = F,
         show_colnames = F, show_rownames = F,
         clustering_method = "ward.D2",
         main = gene, annotation_row = annot.frame,
         annotation_colors = colors_list)
pheatmap((mat.numbers2.dcast2[,c(4,5,6,8)]), scale = "column", #cluster on features
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         breaks = c(-4,seq(-3,3,length=100),4),
         cluster_cols = F, cluster_rows = T,
         show_colnames = T, show_rownames = F,
         clustering_method = "ward.D2",
         main = gene, annotation_row = annot.frame[, 1,drop=F],
         annotation_colors = colors_list)


#multiple genes-----
gene = "Tnf"#"Ifit3"#"Tnf" #M0:"Mx2", "Nfkbiz", "Egr3"
X0.mat.numbers2.dcast = dcast(mat.numbers2, path_stim~time, value.var = gene)
X0.mat.numbers2.dcast$path_stim = gsub('[0-9]+_', '', X0.mat.numbers2.dcast$path_stim)
# X0.mat.numbers2.dcast$path_stim = gsub(".*_","",X0.mat.numbers2.dcast$path_stim)

gene =  "Il6" #"Cmpk2"#"Mx2","Nfkbiz"
X1.mat.numbers2.dcast = dcast(mat.numbers2, path_stim~time, value.var = gene)
X1.mat.numbers2.dcast$path_stim = gsub('[0-9]+_', '', X0.mat.numbers2.dcast$path_stim)

gene = "Ccl5"  #"Irf1"
X2.mat.numbers2.dcast = dcast(mat.numbers2, path_stim~time, value.var = gene)
X2.mat.numbers2.dcast$path_stim = gsub('[0-9]+_', '', X0.mat.numbers2.dcast$path_stim)

gene = "Cxcl10"
X3.mat.numbers2.dcast = dcast(mat.numbers2, path_stim~time, value.var = gene)
X3.mat.numbers2.dcast$path_stim = gsub('[0-9]+_', '', X0.mat.numbers2.dcast$path_stim)

gene = "Nfkbia"#"Gna15","Egr3"
X4.mat.numbers2.dcast = dcast(mat.numbers2, path_stim~time, value.var = gene)
X4.mat.numbers2.dcast$path_stim = gsub('[0-9]+_', '', X4.mat.numbers2.dcast$path_stim)
# X4.mat.numbers2.dcast$path_stim = gsub(".*_","",X4.mat.numbers2.dcast$path_stim)

mat.numbers2.dcast = cbind(X0.mat.numbers2.dcast,
                           X1.mat.numbers2.dcast[,-1], X2.mat.numbers2.dcast[,-1],X3.mat.numbers2.dcast[,-1],
                           X4.mat.numbers2.dcast[,-1])
colors_list = list(path_stim = c(CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3"))
# colors_list = list(path_stim = c(`0_rep2only_TNF`="darkred", `1_IFNg_TNF`="#00BA38", `2_IL4_gt80_TNF`="#619CFF"))
annot.frame = mat.numbers2.dcast[,c(1),drop=F] #c(13,4,5,6,8,9,10)
annot.frame$path_stim = gsub("_..*","", annot.frame$path_stim)

mat.numbers2.dcast2 = mat.numbers2.dcast[grepl("rep2", mat.numbers2.dcast$path_stim),]
mat.numbers2.dcast2$path_stim = gsub("_..*","", mat.numbers2.dcast2$path_stim)
table(mat.numbers2.dcast2$path_stim)
gc();gc();gc()
mat.numbers2.dcast2=mat.numbers2.dcast2[ , colSums(is.na(mat.numbers2.dcast2))==0]
pheatmap(as.matrix(mat.numbers2.dcast2[order(mat.numbers2.dcast2$path_stim),-c(1)]), 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
         breaks = c(0,seq(0.01,0.99,length=100),1),
         cluster_cols = F, cluster_rows = T, # cluster_rows=as.hclust(row_dend),
         show_colnames = F, show_rownames = F,
         annotation_colors = colors_list,
         clustering_method = "ward.D2", cutree_rows = 6, gaps_col = c(102),
         main = gene, annotation_row = annot.frame)


###############################################################################
# Figure 3-4: calculate trajectory features----
# calculate features on single cell RNA trajectories,scaled-----
if(1){ 
  #Calculate peak induction of each gene
  dynamics = data.frame()
  for (i in 1:dim(A)[3]){
    print(i)
    # gene_name = colnames(reconstructed_pc$reconstructed_trajectories)[i]
    gene_name = colnames(mat.numbers)[i]
    print(gene_name)
    A.subset = as.data.frame(t(A[ , ,i]@data))
    my.dataframe = cbind(label = labels, A.subset)
    peak_amp <- apply(my.dataframe[,-1], 1, max) 
    tmp = data.frame(peak_amp =peak_amp, stimulus =labels,type=labels2, gene = gene_name)
    dynamics <- rbind(dynamics, tmp)
  }
  # write.table(dynamics, "./trajectory/trajectory_features_M0rep2only_peakamp_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M1_IFNg_peakamp_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M2_IL4_peakamp_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_allM0M1M2_peakamp_k20.txt", quote=F,row.names = F, sep = "\t")
  
  dynamics = read.delim("./trajectory/trajectory_features_allM0M1M2_peakamp_k20.txt")
  ggplot(dynamics[grepl("Tnf$",dynamics$gene),], aes(stimulus, peak_amp))+
    facet_grid(~type)+
    geom_point(position = "jitter",alpha = 0.5)+geom_violin(aes(color = stimulus), outlier.shape = NA)+ylim(0,1)+
    stat_summary(fun.y = median, geom='point', size = 2, colour = "blue")+
    theme_bw(base_size = 16)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position = "none")
  
  #peak fold change
  dynamics = data.frame()
  for (i in 1:dim(A)[3]){
    print(i)
    # gene_name = colnames(reconstructed_pc$reconstructed_trajectories)[i]
    gene_name = colnames(mat.numbers)[i]
    print(gene_name)
    A.subset = as.data.frame(t(A[ , ,i]@data))
    my.dataframe = cbind(label = labels, A.subset)
    peak_amp <- apply(my.dataframe[,-1], 1, max) 
    tmp = data.frame(peak_amp_lfc = log2((peak_amp/(my.dataframe$V1+0.01))+1), 
                     peak_amp_fc = (peak_amp/(my.dataframe$V1+0.01)),
                     time0_amp = my.dataframe$V1,
                     stimulus =labels, type=labels2,gene = gene_name)
    dynamics <- rbind(dynamics, tmp)
  }
  # write.table(dynamics, "./trajectory/trajectory_features_M0rep2only_peakamplogFC_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M1_IFNg_peakamplogFC_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M2_IL4_peakamplogFC_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_allM0M1M2_peakamplogFC_k20.txt", quote=F,row.names = F, sep = "\t")
  
  # dynamics = read.delim("./trajectory/trajectory_features_allM0M1M2_peakamplogFC_k20.txt")
  ggplot(dynamics[grepl("Tnf$",dynamics$gene),], aes(stimulus, peak_amp_lfc))+
    geom_point(position = "jitter", alpha = 0.5)+geom_violin(aes(color = stimulus), outlier.shape = NA)+
    stat_summary(fun.y = median, geom='point', size = 2, colour = "blue")+
    theme_bw(base_size = 16)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
  
  
  #Speed at time 1hr
  timept_tangent = num_timepts/8 +2 #for 1hr
  dynamics = data.frame()
  for (i in 1:dim(A)[3]){
    print(i)
    # gene_name = colnames(reconstructed_pc$reconstructed_trajectories)[i]
    gene_name = colnames(mat.numbers)[i]
    print(gene_name)
    A.subset = as.data.frame(t(A[ , ,i]@data))
    my.dataframe = cbind(label = labels, A.subset)
    
    colnames(my.dataframe)[-1] <- unique(mat.meta$time)
    timeseg <- as.numeric(names(my.dataframe[,-1])[round(timept_tangent)+2]) - 
      as.numeric(names(my.dataframe[,-1])[round(timept_tangent)-2])
    rise <- my.dataframe[,round(timept_tangent)+2]- my.dataframe[,round(timept_tangent)-2]
    
    tmp = data.frame(speed1hr = (rise/timeseg), stimulus =labels, type=labels2,gene = gene_name)
    dynamics <- rbind(dynamics, tmp)
  }
  # write.table(dynamics, "./trajectory/trajectory_features_M0rep2only_speed1hr_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M1_IFNg_speed1hr_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M2_IL4_speed1hr_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_allM0M1M2_speed1hr_k20.txt", quote=F,row.names = F, sep = "\t")
  
  ggplot(dynamics[grepl("Tnf$",dynamics$gene),], aes(stimulus, speed1hr))+
    geom_point(position = "jitter",alpha=0.4)+geom_violin(aes(color = stimulus), outlier.shape = NA)+
    stat_summary(fun.y = median, geom='point', size = 2, colour = "blue")+
    theme_bw(base_size = 16)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
  
  
  # Integral, total mRNA
  dynamics = data.frame()
  for (i in 1:dim(A)[3]){
    tryCatch(
      expr = {
        print(i)
        # gene_name = colnames(reconstructed_pc$reconstructed_trajectories)[i]
        gene_name = colnames(mat.numbers)[i]
        print(gene_name)
        A.subset = as.data.frame(t(A[ , ,i]@data))
        my.dataframe = cbind(label = labels, A.subset)
        colnames(my.dataframe)[-1] <- unique(mat.meta$time)
        
        time <- unique(mat.meta$time)
        integral <- apply(my.dataframe[,-1], 1, function(x) unlist(integrate(approxfun(time, x), range(time)[1], range(time)[2],rel.tol =.Machine$double.eps^.2))$value)
        
        tmp = data.frame(integral = integral, stimulus =labels, type=labels2,gene = gene_name)
        dynamics <- rbind(dynamics, tmp)
      }, error = function(e){
        message('Caught an error!')
        integral <- NA
        tmp = data.frame(integral = integral, stimulus =labels, type=labels2,gene = gene_name)
        dynamics <- rbind(dynamics, tmp)
      }
    )
  }
  # write.table(dynamics, "./trajectory/trajectory_features_M0rep2only_integral_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M1_IFNg_integral_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M2_IL4_integral_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_allM0M1M2_integral_k20.txt", quote=F,row.names = F, sep = "\t")
  
  # dynamics = read.delim("./trajectory/trajectory_features_allM0M1M2_integral_k20.txt")
  ggplot(dynamics[grepl("Nfkbia$",dynamics$gene),], aes(stimulus, integral))+
    geom_violin(aes(color = stimulus), outlier.shape = NA)+geom_point(position = "jitter", alpha = 0.5, size=0.1)+
    facet_grid(~type)+
    stat_summary(fun.y = median, geom='point', size = 2, colour = "blue")+
    theme_bw(base_size = 16)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
  
}
if(0){
  #Speed to peak
  dynamics = data.frame()
  for (i in 1:dim(A)[3]){
    print(i)
    gene_name = colnames(reconstructed_pc$reconstructed_trajectories)[i]
    print(gene_name)
    A.subset = as.data.frame(t(A[ , ,i]@data))
    my.dataframe = cbind(label = labels, A.subset)
    
    colnames(my.dataframe)[-1] <- unique(mat.meta$time)
    time2max <- apply(my.dataframe[,-1], 1, function(x) as.numeric(names(x)[which(x == max(x))]))
    logFCmax <- apply(my.dataframe[,-1], 1, max) - my.dataframe[,-1][,1]
    
    tmp = data.frame(speed = (logFCmax/time2max), time2peak = time2max, logFCmax = logFCmax, stimulus =labels, gene = gene_name)
    dynamics <- rbind(dynamics, tmp)
  }
  # write.table(dynamics, "./trajectory/trajectory_features_M0rep2only_speed2peak_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M1_IFNg_speed2peak_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M2_IL4_speed2peak_k20.txt", quote=F,row.names = F, sep = "\t")
  ggplot(dynamics[grepl("Tnf$",dynamics$gene),], aes(stimulus, speed))+
    geom_point(position = "jitter",alpha=0.4)+geom_violin(aes(color = stimulus), outlier.shape = NA)+
    stat_summary(fun.y = median, geom='point', size = 2, colour = "blue")+
    theme_bw(base_size = 16)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
  ggplot(dynamics[grepl("Tnf$",dynamics$gene),], aes(stimulus, time2peak))+
    geom_point(position = "jitter",alpha=0.4)+geom_violin(aes(color = stimulus), outlier.shape = NA)+
    stat_summary(fun.y = median, geom='point', size = 2, colour = "blue")+
    theme_bw(base_size = 16)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
  
  
  #duration >0.2
  8/num_timepts*60 #~4.6minutes per frame
  dynamics = data.frame()
  for (i in 1:dim(A)[3]){
    print(i)
    gene_name = colnames(reconstructed_pc$reconstructed_trajectories)[i]
    print(gene_name)
    A.subset = as.data.frame(t(A[ , ,i]@data))
    my.dataframe = cbind(label = labels, A.subset)
    
    colnames(my.dataframe)[-1] <- unique(mat.meta$time)
    above0.1 <- apply(my.dataframe[,-1], 1, function(x){sum((x>0.2)==T)} )
    
    tmp = data.frame(duration = (8/num_timepts) *above0.1, stimulus =labels, gene = gene_name)
    dynamics <- rbind(dynamics, tmp)
  }
  # write.table(dynamics, "./trajectory/trajectory_features_M0rep2only_duration_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M1_IFNg_duration_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M2_IL4_duration_k20.txt", quote=F,row.names = F, sep = "\t")
  ggplot(dynamics[grepl("Tnf$",dynamics$gene),], aes(stimulus, duration))+
    geom_point(position = "jitter",alpha=0.4)+geom_violin(aes(color = stimulus), outlier.shape = NA)+
    stat_summary(fun.y = median, geom='point', size = 2, colour = "blue")+
    theme_bw(base_size = 16)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
  
  #duration >0.3, or 0.25
  8/num_timepts*60 #~4.6minutes per frame
  dynamics = data.frame()
  for (i in 1:dim(A)[3]){
    print(i)
    gene_name = colnames(reconstructed_pc$reconstructed_trajectories)[i]
    print(gene_name)
    A.subset = as.data.frame(t(A[ , ,i]@data))
    my.dataframe = cbind(label = labels, A.subset)
    
    colnames(my.dataframe)[-1] <- unique(mat.meta$time)
    above0.1 <- apply(my.dataframe[,-1], 1, function(x){sum((x>0.25)==T)} )
    
    tmp = data.frame(duration = (8/num_timepts) *above0.1, stimulus =labels, gene = gene_name)
    dynamics <- rbind(dynamics, tmp)
  }
  # write.table(dynamics, "./trajectory/trajectory_features_M0rep2only_duration0.3_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M1_IFNg_duration0.3_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M2_IL4_duration0.3_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M0rep2only_duration0.25_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M1_IFNg_duration0.25_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M2_IL4_duration0.25_k20.txt", quote=F,row.names = F, sep = "\t")
  ggplot(dynamics[grepl("Tnf$",dynamics$gene),], aes(stimulus, duration))+
    geom_point(position = "jitter",alpha=0.4)+geom_violin(aes(color = stimulus), outlier.shape = NA)+
    stat_summary(fun.y = median, geom='point', size = 2, colour = "blue")+
    theme_bw(base_size = 16)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
  
  
  #time2halfintegral
  dynamics = data.frame()
  for (i in 1:dim(A)[3]){
    tryCatch(
      expr = {
        print(i)
        gene_name = colnames(reconstructed_pc$reconstructed_trajectories)[i]
        print(gene_name)
        A.subset = as.data.frame(t(A[ , ,i]@data))
        my.dataframe = cbind(label = labels, A.subset)
        
        colnames(my.dataframe)[-1] <- unique(mat.meta$time)
        time <- unique(mat.meta$time)
        
        tmp <- apply(my.dataframe[-1], 1, function(x) {
          halfmax <- unlist(integrate(approxfun(time, x), range(time)[1], range(time)[2], rel.tol =.Machine$double.eps^.2))$value / 2
          root <- function(y) {
            return(unlist(integrate(approxfun(time, x), 0, y,rel.tol =.Machine$double.eps^.2))$value - halfmax)
          }
          return(uniroot(root, c(0, 8))$root)
        })
        
        tmp = data.frame(earlyVSlate = tmp, stimulus =labels, gene = gene_name)
        dynamics <- rbind(dynamics, tmp)
      }, error = function(e){
        message('Caught an error!')
        tmp <- NA
        tmp = data.frame(earlyVSlate = tmp, stimulus =labels, gene = gene_name)
        dynamics <- rbind(dynamics, tmp)
      }
    )
  }
  # write.table(dynamics, "./trajectory/trajectory_features_M0rep2only_time2halfint_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M1_IFNg_time2halfint_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M2_IL4_time2halfint_k20.txt", quote=F,row.names = F, sep = "\t")
  
  #peak repression FC
  dynamics = data.frame()
  for (i in 1:dim(A)[3]){
    print(i)
    gene_name = colnames(reconstructed_pc$reconstructed_trajectories)[i]
    print(gene_name)
    A.subset = as.data.frame(t(A[ , ,i]@data))
    my.dataframe = cbind(label = labels, A.subset)
    peak_amp <- apply(my.dataframe[,-1], 1, max) 
    
    colnames(my.dataframe)[-1] <- unique(mat.meta$time)
    my.dataframe$time2max <- apply(my.dataframe[,-1], 1, function(x) (names(x)[which(x == max(x))]))
    trough <- apply(my.dataframe[,-1], 1, function(x) 
      min(as.numeric(x[ which(names(x) == x["time2max"]) : which(names(x) == "8") ])))#trough after peak time
    trough_alltime = apply(my.dataframe[,-1], 1, min) #min across all time
    tmp = data.frame(trough_postpeak = trough,
                     trough_alltime = trough_alltime,
                     peak_repress_lfc = log2((peak_amp/trough)+1), 
                     peak_repress_fc = (peak_amp/trough),
                     stimulus =labels, gene = gene_name)
    dynamics <- rbind(dynamics, tmp)
  }
  # write.table(dynamics, "./trajectory/trajectory_features_M0rep2only_peakrepresslogFC_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M1_IFNg_peakrepresslogFC_k20.txt", quote=F,row.names = F, sep = "\t")
  # write.table(dynamics, "./trajectory/trajectory_features_M2_IL4_peakrepresslogFC_k20.txt", quote=F,row.names = F, sep = "\t")
  dynamics = read.delim("./trajectory/trajectory_features_M0rep2only_peakrepresslogFC_k20.txt")
  ggplot(dynamics[grepl("Nfkbia$",dynamics$gene),], aes(stimulus, peak_repress_lfc))+
    geom_point(position = "jitter", alpha = 0.5)+geom_violin(aes(color = stimulus), outlier.shape = NA)+
    stat_summary(fun.y = median, geom='point', size = 2, colour = "blue")+
    theme_bw(base_size = 16)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")
  
  
  
}
##############################################################################
# read in dynamical features-----
#read data
dynamics1 = read.delim("./trajectory/trajectory_features_M0rep2only_peakamp_k20.txt");dynamics1$cell=rep(1:5996)
dynamics2 = read.delim("./trajectory/trajectory_features_M0rep2only_peakamplogFC_k20.txt")
dynamics3 = read.delim("./trajectory/trajectory_features_M0rep2only_duration0.25_k20.txt")
dynamics4 = read.delim("./trajectory/trajectory_features_M0rep2only_integral_k20.txt")
dynamics5 = read.delim("./trajectory/trajectory_features_M0rep2only_speed2peak_k20.txt")
dynamics6 = read.delim("./trajectory/trajectory_features_M0rep2only_speed1hr_k20.txt")
dynamics7 = read.delim("./trajectory/trajectory_features_M0rep2only_time2halfint_k20.txt");dynamics7$cell=rep(1:5996)
dynamics8 = read.delim("./trajectory/trajectory_features_M0rep2only_peakrepresslogFC_k20.txt")
dynamics_M0 = cbind(cell = dynamics1$cell, dynamics1[, c(3,2,1)], peak_amp_lfc=dynamics2$peak_amp_lfc,integral=dynamics4$integral,
                    time2peak = dynamics5$time2peak,
                    speed1hr=dynamics6$speed1hr,duration0.25=dynamics3$duration, 
                    time2halfint=dynamics7$earlyVSlate, 
                    trough_postpeak = dynamics8$trough_postpeak,
                    peak_repress_lfc=dynamics8$peak_repress_lfc)
dynamics_M0$peak_repress_lfc = log2(dynamics_M0$peak_amp/(dynamics_M0$trough_postpeak+0.01)+1)

dynamics1 = read.delim("./trajectory/trajectory_features_M1_IFNg_peakamp_k20.txt");dynamics1$cell=rep(1:6000)
dynamics2 = read.delim("./trajectory/trajectory_features_M1_IFNg_peakamplogFC_k20.txt")
dynamics3 = read.delim("./trajectory/trajectory_features_M1_IFNg_duration0.25_k20.txt")
dynamics4 = read.delim("./trajectory/trajectory_features_M1_IFNg_integral_k20.txt");dynamics4$cell=rep(1:6000)
dynamics5 = read.delim("./trajectory/trajectory_features_M1_IFNg_speed2peak_k20.txt")
dynamics6 = read.delim("./trajectory/trajectory_features_M1_IFNg_speed1hr_k20.txt")
dynamics7 = read.delim("./trajectory/trajectory_features_M1_IFNg_time2halfint_k20.txt");dynamics7$cell=rep(1:6000)
dynamics8 = read.delim("./trajectory/trajectory_features_M1_IFNg_peakrepresslogFC_k20.txt")
dynamics_M1 = cbind(cell = dynamics1$cell, dynamics1[, c(3,2,1)], peak_amp_lfc=dynamics2$peak_amp_lfc,
                    integral=dynamics4$integral,
                    time2peak = dynamics5$time2peak,
                    speed1hr=dynamics6$speed1hr,duration0.25=dynamics3$duration,
                    time2halfint=dynamics7$earlyVSlate, 
                    trough_postpeak = dynamics8$trough_postpeak,
                    peak_repress_lfc=dynamics8$peak_repress_lfc)
dynamics_M1$peak_repress_lfc = log2(dynamics_M1$peak_amp/(dynamics_M1$trough_postpeak+0.01)+1)


dynamics1 = read.delim("./trajectory/trajectory_features_M2_IL4_peakamp_k20.txt");dynamics1$cell=rep(1:6000)
dynamics2 = read.delim("./trajectory/trajectory_features_M2_IL4_peakamplogFC_k20.txt")
dynamics3 = read.delim("./trajectory/trajectory_features_M2_IL4_duration0.25_k20.txt")
dynamics4 = read.delim("./trajectory/trajectory_features_M2_IL4_integral_k20.txt");dynamics4$cell=rep(1:6000)
dynamics5 = read.delim("./trajectory/trajectory_features_M2_IL4_speed2peak_k20.txt")
dynamics6 = read.delim("./trajectory/trajectory_features_M2_IL4_speed1hr_k20.txt")
dynamics7 = read.delim("./trajectory/trajectory_features_M2_IL4_time2halfint_k20.txt");dynamics7$cell=rep(1:6000)
dynamics8 = read.delim("./trajectory/trajectory_features_M2_IL4_peakrepresslogFC_k20.txt")
dynamics_M2 = cbind(cell = dynamics1$cell, dynamics1[, c(3,2,1)], peak_amp_lfc=dynamics2$peak_amp_lfc,
                    integral=dynamics4$integral,
                    time2peak = dynamics5$time2peak,
                    speed1hr=dynamics6$speed1hr,duration0.25=dynamics3$duration,
                    time2halfint=dynamics7$earlyVSlate,#[match(paste0(dynamics1$gene, dynamics1$stimulus, dynamics1$cell),paste0(dynamics7$gene, dynamics7$stimulus, dynamics7$cell))], 
                    trough_postpeak = dynamics8$trough_postpeak,
                    peak_repress_lfc=dynamics8$peak_repress_lfc)
dynamics_M2$peak_repress_lfc = log2(dynamics_M2$peak_amp/(dynamics_M2$trough_postpeak+0.01)+1)

#################### read data scaled across all macstates----------------
dynamics1 = read.delim("./trajectory/trajectory_features_allM0M1M2_peakamp_k20.txt");dynamics1$cell=rep(1:17996)
dynamics2 = read.delim("./trajectory/trajectory_features_allM0M1M2_peakamplogFC_k20.txt")
dynamics3 = read.delim("./trajectory/trajectory_features_allM0M1M2_speed1hr_k20.txt")
dynamics4 = read.delim("./trajectory/trajectory_features_allM0M1M2_integral_k20.txt")
dynamics_M0M1M2 = cbind(cell = dynamics1$cell, dynamics1[, c(4,3,2,1)], 
                        peak_amp_lfc=dynamics2$peak_amp_lfc,
                        speed1hr=dynamics3$speed1hr,
                        integral=dynamics4$integral)

###############################################################################
# Figure 3-4: calculate mutual information on trajectory features----
# mutual info on dynamical features for each gene (subtract dynamic vs static)----
library(SLEMI)
dynamics = dynamics_M0
dynamics = dynamics_M1
dynamics = dynamics_M2
head(dynamics)
collect_mi = data.frame()

# for (i in c( "integral", "duration0.25", "time2halfint","peak_repress_lfc", "time2peak","duration")){
# for (i in c("integral","peak_amp","peak_amp_lfc", "speed1hr")){
  for (g in unique(dynamics$gene)){
    print(i)
    print(g)
    my.dataframe.subset = dynamics[dynamics$gene==g, c("gene","stimulus", i)]
    
    tryCatch(
      expr = {
        output_capacity <- capacity_logreg_main(my.dataframe.subset, signal = "stimulus", response = i,
                                                output_path = NULL, testing = F)
        
        tmp = data.frame(feature = i,  gene = g, cc = output_capacity$cc)
        collect_mi = rbind(collect_mi, tmp)
      }, error = function(e){
        message("Error for gene!")
        tmp = data.frame(feature = i,  gene = g, cc = NA)
        collect_mi = rbind(collect_mi, tmp)
      }
      
    )
  
}
# write.table(collect_mi, "./infotheo/dynamics_k20_trajectoryfeatures_singlegene_SLEMI_M0.txt", row.names = F, sep = "\t", quote=F)
# write.table(collect_mi, "./infotheo/dynamics_k20_trajectoryfeatures_singlegene_SLEMI_M1.txt", row.names = F, sep = "\t", quote=F)
# write.table(collect_mi, "./infotheo/dynamics_k20_trajectoryfeatures_singlegene_SLEMI_M2.txt", row.names = F, sep = "\t", quote=F)

#####################################################################################
# Figure 4: collect dim best start with 1D top 20 (done on Precision)---- 
collect_all = read.delim("./infotheo/dynamics_k20_trajectoryfeatures_singlegene_SLEMI_M0.txt")
collect_all = read.delim("./infotheo/dynamics_k20_trajectoryfeatures_singlegene_SLEMI_M1.txt")
collect_all = read.delim("./infotheo/dynamics_k20_trajectoryfeatures_singlegene_SLEMI_M2.txt")

dynamics = dynamics_M0
dynamics = dynamics_M1
dynamics = dynamics_M2

collect_all = collect_all[order(collect_all$cc, decreasing = T), ]
require(data.table) ## 1.9.2
collect_all <- as.data.table(collect_all)
collect_all_best = collect_all[collect_all[, .I[cc == max(cc)], by="feature"]$V1] 
collect_all_best$dim = 0
collect_all_best = collect_all_best[, c("feature", "dim","cc",  "gene")]
table(collect_all_best$feature)

collect = data.frame()
collect_dimensionbest = data.frame()

for (i in c("peak_amp","peak_amp_lfc", "integral", "speed", "time2peak", "logFCmax", "speed1hr", "duration")){
  print(i)
  
  collect_dimension = (collect_all[grepl(i, collect_all$feature),][(1:20),]) #start 1D
  
  my.dataframe = dynamics[, c("gene","stimulus", i)]
  my.dataframe = cbind(cell = rep(1: 5996), my.dataframe)
  my.dataframe = dcast(my.dataframe, stimulus+cell~gene, value.var = i)
  
  for (d in seq(1:5)){
    print(paste0("dimension: ",d))
    
    collect_dimension = collect_dimension[order(collect_dimension$cc, decreasing =T),]
    genesets = collect_dimension$gene[c(1:20)]
    print(genesets)
    
    collect_dimension = data.frame() #start over once got the top20
    for (g in 1:length(genesets)){
      genes = genesets[[g]]
      print(genes)
      
      other_genes = c(colnames(my.dataframe)[!colnames(my.dataframe) %in% genes])
      other_genes
      for (a in 1:length(other_genes)){
        added_gene = other_genes[a+1]
        added_gene
        my.dataframe.subset = my.dataframe[, colnames(my.dataframe) %in% c("stimulus", genes, added_gene)]
        str(my.dataframe.subset)
        
        output_capacity <- capacity_logreg_main(my.dataframe.subset, signal = "stimulus", response = colnames(my.dataframe.subset)[-1],
                                                output_path = NULL, testing = F)
        
        tmp = data.frame(feature = i,  dim = d, cc = output_capacity$cc)
        tmp$gene = list(c(genes, added_gene))
        collect_dimension = rbind(collect_dimension, tmp)
      }
      
    }
    collect_dimension = collect_dimension[order(collect_dimension$cc, decreasing =T),]
    collect_dimensionbest = rbind(collect_dimensionbest,
                                  data.frame(collect_dimension[1,]))
    
  }
}

# saveRDS(collect_dimensionbest,"./infotheo/collect_dimensionbest_dynamics_trajfeatures_M0_May2022.rds")

#######################################################################
# Figure 5: gene-gene correlations----
gene = "Nfkbia$|Tnf$"
dynamics = rbind(cbind(dynamics_M0, type = "M0"),
                 cbind(dynamics_M1, type = "M1"),
                 cbind(dynamics_M2, type = "M2")) 
# dynamics = dynamics_M0M1M2
dynamics.dcast_lfc= dcast(dynamics, type+cell+stimulus ~ gene, value.var = "peak_amp_lfc")
dynamics.dcast_int= dcast(dynamics, type+cell+stimulus ~ gene, value.var = "integral")
dynamics.dcast_amp= dcast(dynamics, type+cell+stimulus ~ gene, value.var = "peak_amp")
dynamics.dcast_spd= dcast(dynamics, type+cell+stimulus ~ gene, value.var = "speed1hr")
dynamics.dcast = cbind(lfc =dynamics.dcast_lfc, int=dynamics.dcast_int[,-c(1:3)],
                       amp =dynamics.dcast_amp[,-c(1:3)], spd=dynamics.dcast_spd[,-c(1:3)])
colnames(dynamics.dcast)[1:3] = gsub("lfc.","", colnames(dynamics.dcast)[1:3])
colnames(dynamics.dcast)
ggplot(dynamics.dcast[grepl("",dynamics.dcast$type),], aes(int.Tnf, int.Cxcl10 ))+
  facet_wrap(~type, scales = "free",nrow = 1)+
  stat_cor(method = "pearson", label.x = 1)+
  # ylim(1,3.5)+ #for lfc plotting
  # stat_summary(aes(amp.Tnf, lfc.Tnf, group = stimulus,colour=stimulus), fun.y=mean,  geom="point")+
  geom_smooth( method = "lm", alpha=0.1)+
  geom_point(aes(color = stimulus),alpha=0.3)+theme_bw(base_size = 18)+ggtitle(gene)+theme(legend.position = "none")
ggplot(dynamics.dcast[grepl("CpG|LPS|TNF",dynamics.dcast$stimulus),], aes(spd.Tnf, spd.Cxcl10 ))+
  facet_wrap(~stimulus, scales = "free",nrow = 2)+
  geom_smooth(aes(group=type, color = type), method = "lm", alpha=0.1)+
  geom_point(aes(color = type),alpha=0.5, size=0.5)+
  stat_cor(aes(color =type),method = "pearson", label.x = 0)+
  theme_bw(base_size = 14)+ggtitle(gene)+theme(legend.position = "none")

#make correlation violin/bar plots------
combine = data.frame()
feature = "spd"
for (t in c("M0","M1","M2")){
  for (s in c("CpG","IFNb","LPS","P3CSK","PIC","TNF")){
    print(t);print(s)
    corr = cor(dynamics.dcast[grepl(t,dynamics.dcast$type)&grepl(s,dynamics.dcast$stimulus),grepl(feature,colnames(dynamics.dcast))])
    tmp = melt(corr)
    tmp$type = t
    tmp$stimulus = s
    combine = rbind(combine, tmp)
  }
}
head(combine)
combine$Var1 = gsub(paste0(feature,"."),"", combine$Var1)
combine$Var2 = gsub(paste0(feature,"."),"", combine$Var2)
clusters = readxl::read_excel("F://scRNAseq_macro/SuppTables/TableS4_gene_regulatory_strategies_allgenes.xlsx")
combine$GRS1 = clusters$clusters[match(combine$Var1, clusters$gene)]
combine$GRS2 = clusters$clusters[match(combine$Var2, clusters$gene)]
combine$sameGRS = ifelse(combine$GRS1==combine$GRS2, "same", "diff")

combine$value = ifelse(combine$value==1,NA, combine$value)
geneset = "Tnf$|Il1b"
tmp =na.omit(combine[grepl(geneset,combine$Var1)&grepl(geneset, combine$Var2)&grepl("", combine$GRS1),])
ggplot(na.omit((combine[grepl(geneset,combine$Var1)&grepl(geneset, combine$Var2)&grepl("", combine$GRS1),])), 
       aes(stimulus, value/2))+
  geom_bar(stat = "identity", fill = "palegreen")+
  facet_grid(~type)+ ylab("Correlation Coeff.")+ggtitle(geneset)+ylim(0,.75)+
  # geom_boxplot(outlier.shape = NA)+geom_point(aes(color = GRS1),size=0.01,position = "jitter")+
  theme_bw(base_size = 14)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


#######################################################################
# Figure 6: response specificity profile----
#PCA on dynamical features, M0M1M2 combined----
dynamics_M0$cell = paste0(dynamics_M0$cell, "_M0")
dynamics_M1$cell = paste0(dynamics_M1$cell, "_M1")
dynamics_M2$cell = paste0(dynamics_M2$cell, "_M2")
dynamics = rbind(dynamics_M0, dynamics_M1, dynamics_M2)
dynamics$type = gsub(".*_", "", dynamics$cell)
for (feature in c("peak_amp","peak_amp_lfc" ,"integral" ,"speed1hr","duration0.25","time2halfint")){
  dynamics.dcast = dcast(dynamics, gene~cell, value.var = feature)
  write.table(na.omit(dynamics.dcast), paste0("./trajectory/matrix_trajectory_features_M0M1M2_",feature,"_k20.txt"),
              row.names = F, sep = "\t", quote=F)
}
for (feature in c("peak_amp","peak_amp_lfc" ,"integral" ,"speed1hr","duration0.25","time2halfint")){
  PCA_from_file(paste0("./trajectory/matrix_trajectory_features_M0M1M2_",feature,"_k20.txt"),
                center =T, scale = F)
}
for (feature in c("peak_amp","peak_amp_lfc" ,"integral" ,"speed1hr","duration0.25","time2halfint")){
  setwd("F://scRNAseq_macro/scRNAseq_macro/trajectory/")
  plot_pca(paste0("./matrix_trajectory_features_M0M1M2_",feature,"_k20_prcomp_scores.txt"),
           info.name = paste0("X",dynamics$cell), info.type = as.factor(dynamics$stimulus), labels = F,
           alpha=0.5, pt.size=.5, title = feature)
  plot_pca(paste0("./matrix_trajectory_features_M0M1M2_",feature,"_k20_prcomp_scores.txt"),
           info.name = paste0("X",dynamics$cell), info.type = as.factor(dynamics$type), labels = F,
           alpha=0.5, pt.size=.5, title = feature)
  setwd("F://scRNAseq_macro/scRNAseq_macro/")
}
# plot PCA loadings-------
feature = "peak_amp_lfc"#"speed1hr" #"peak_amp_lfc"
loadings = read.delim(paste0("./trajectory/matrix_trajectory_features_M0M1M2_",feature,"_k20_prcomp_loadings.txt"))
clusters = readxl::read_excel("F://scRNAseq_macro/SuppTables/TableS4_gene_regulatory_strategies_allgenes.xlsx")
loadings$GRS = clusters$clusters[match(loadings$Loading, clusters$gene)]
ggplot(loadings, aes(PC1, PC2))+geom_point(aes(color = GRS), size=2)+
  theme_bw(base_size = 16)+ggtitle(feature)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)

#pairwise max MI on PC components, all M0M1M2 together----
collect = data.frame()
nPCs = 3
list = c("CpG", "IFNb", "LPS", "P3CSK","PIC", "TNF") 
for (feature in c("peak_amp","peak_amp_lfc" ,"integral" ,"speed1hr","duration0.25","time2halfint")){
  print(feature)
  pc.scores = read.delim(paste0("./trajectory/matrix_trajectory_features_M0M1M2_",feature,"_k20_prcomp_scores.txt"))
  head(pc.scores)[1:5]
  colnames(pc.scores)[1]="Sample"
  pc.scores$stimulus = dynamics$stimulus[match(pc.scores$Sample, paste0("X",dynamics$cell))]
  pc.scores$type = gsub(".*_", "", pc.scores$Sample)
  colnames(pc.scores)
  
  mi.frame = rbind(pc.scores)
  rownames(mi.frame)= mi.frame$Sample
  mi.frame$Sample = mi.frame$stimulus
  colnames(mi.frame)[1]="label"
  
  for (m in c("M0","M1","M2")){
    print(m)
    
    library(fpc);library(philentropy)
    # pca.macro = as.data.frame(Embeddings(macro[["pca"]]))
    
    for (i in seq(1:(length(list)-1) )){
      for (j in seq(i+1,length(list)) ){
        
        
        print(list[i])
        print(list[j])
        
        wanted.1 = rownames(mi.frame)[mi.frame$type == m & mi.frame$stimulus == list[i] ]
        wanted.2 = rownames(mi.frame)[mi.frame$type == m & mi.frame$stimulus == list[j] ]
        
        my.dataframe= mi.frame[(rownames(mi.frame)%in%wanted.1 | rownames(mi.frame)%in%wanted.2), c(1:(nPCs+1))]
        
        # my.dataframe=my.dataframe[grepl(s, my.dataframe$label),]
        #--------------------------------mi using SLEMI -----------------------
        
        #calc information capacity: capacity_logreg_main(dataRaw, signal, response, output_path)
        output_capacity <- capacity_logreg_main(my.dataframe, signal = "label", response = colnames(my.dataframe)[-1],
                                                # paste0("F:/scRNAseq_macro/scRNAseq_macro/infotheo/cc_", i), 
                                                testing = T, boot_num = 10, boot_prob = .5 ,testing_cores = 4)
        sd = sd(sapply(output_capacity$testing$bootstrap, '[[', 3)) #get cc from each bootstrap
        tmp = data.frame(type = m, feature=feature, stim1 = list[i], stim2 = list[j], cc = output_capacity$cc, sd = sd)
        collect = rbind(collect, tmp)
        
      }
    }
  }
}
# write.table(collect, "./infotheo/channel_capacity_pairwise_M0M1M2_reconstruction_integral_RSI_3comps.txt",quote=F, sep="\t",row.names = F)
# write.table(collect, "./infotheo/channel_capacity_pairwise_M0M1M2_reconstruction_integral_RSI_10comps.txt",quote=F, sep="\t",row.names = F)
# write.table(collect, "./infotheo/channel_capacity_pairwise_M0M1M2_reconstruction_integral_RSI_20comps.txt",quote=F, sep="\t",row.names = F)

collect = read.delim("./infotheo/channel_capacity_pairwise_M0M1M2_reconstruction_integral_RSI_3comps.txt")
collect$type = factor(collect$type, levels= c("M0","M1","M2"))
colors_list = (c(M0="#F8766D",M1="#00BA38",M2="#619CFF"))
collect <- transform(collect, stim1 = ifelse(stim1>stim2, stim2, stim1), 
                     stim2 = ifelse(stim1>stim2, stim1,stim2))
collect$pair = paste0(collect$stim1, "_",collect$stim2)
ggplot(na.omit(collect[grepl("", collect$pair),]), aes(pair, cc, fill = type))+
  geom_bar(stat="identity", position ="dodge")+theme_bw(base_size = 12)+ylab("max MI (bits)")+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.5,position=position_dodge(1))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_fill_manual(values=colors_list)+
  facet_wrap(~feature)+
  geom_hline(yintercept = 1, linetype="dotted")+
  # geom_hline(yintercept = log(3)/log(2), linetype="dotted")+
  geom_hline(yintercept = 2, linetype="dotted")+ylim(0,1)

#plot matrix form
ggplot(na.omit(collect[grepl("M0|M1|M2", collect$type),]), aes(pair, cc, fill = type))+
  geom_bar(stat="identity", position ="dodge")+theme_bw(base_size = 12)+ylab("max MI (bits)")+
  geom_errorbar(aes(ymin=cc-sd, ymax=cc+sd), width=.5,position=position_dodge(1))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_fill_manual(values=colors_list)+ 
  facet_grid(type~feature)+
  geom_hline(yintercept = 1, linetype="dotted")+
  # geom_hline(yintercept = log(3)/log(2), linetype="dotted")+
  geom_hline(yintercept = 0.5, linetype="dotted")+ylim(0,1)

#plot subtraction profile----
collect$diff = collect$cc - c(rep(collect$cc[grepl("M0",collect$type)&grepl("peak_amp$", collect$feature)] ,3),
                              rep(collect$cc[grepl("M0",collect$type)&grepl("peak_amp_lfc", collect$feature)] ,3),
                              rep(collect$cc[grepl("M0",collect$type)&grepl("integ", collect$feature)] ,3),
                              rep(collect$cc[grepl("M0",collect$type)&grepl("speed1hr", collect$feature)] ,3),
                              rep(collect$cc[grepl("M0",collect$type)&grepl("durat", collect$feature)] ,3),
                              rep(collect$cc[grepl("M0",collect$type)&grepl("time2halfint", collect$feature)] ,3))
ggplot(na.omit(collect[grepl("peak_amp|integral|speed",collect$feature),]), aes(pair, diff, fill = type))+
  geom_bar(stat="identity", position ="dodge")+theme_bw(base_size = 10)+ylab("max info. difference")+
  geom_errorbar(aes(ymin=diff-sd, ymax=diff+sd), width=.5,position=position_dodge(1))+
  facet_grid(type~feature)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_fill_manual(values=colors_list)

if(1){
  collect = read.delim("./infotheo/channel_capacity_pairwise_M0M1M2_reconstruction_integral_RSI_3comps.txt")
  collect$type = factor(collect$type, levels= c("M0","M1","M2"))
  collect$pair = paste0(collect$stim1, "_",collect$stim2)
  tmp = na.omit(collect)
}
tmp.cast=  dcast(tmp, pair+feature~type, value.var = 'cc')
tmp.cast=  cbind(tmp.cast, sd = (dcast(tmp, pair+feature~type, value.var = 'sd')[,-c(1:2)])^2)# convert to variances

###plot with error bars----
tmp.cast.feature = tmp.cast
tmp.cast.feature$M0.score = (tmp.cast.feature$M0 - tmp.cast.feature$M0)^2
tmp.cast.feature$M1.score = (tmp.cast.feature$M1 - tmp.cast.feature$M0)^2
tmp.cast.feature$M2.score = (tmp.cast.feature$M2 - tmp.cast.feature$M0)^2

tmp.cast.feature$M0.sd = (tmp.cast.feature$sd.M0 + tmp.cast.feature$sd.M0)
tmp.cast.feature$M1.sd = (tmp.cast.feature$sd.M1 + tmp.cast.feature$sd.M0)
tmp.cast.feature$M2.sd = (tmp.cast.feature$sd.M2 + tmp.cast.feature$sd.M0)

tmp.cast.feature$feature = as.character(tmp.cast.feature$feature)
score = data.frame(score = (aggregate(tmp.cast.feature[,c(9:11)], by = list(tmp.cast.feature$feature),FUN = sum) ),
                   sd = (aggregate(tmp.cast.feature[,c(12:14)], by = list(tmp.cast.feature$feature), FUN = sum) ))
score = cbind(feature = score$score.Group.1, sqrt(score[,c(2:4,6:8)]))
score.m = cbind(melt(score[,1:4]), sd=melt(score[,c(1,5:7)])[,3] )
# score$type = factor(score$type, levels= c("M0.score","M1.score","M2.score" ))
# colors_list = (c(M0.score="#F8766D",M1.score="#00BA38",M2.score="#619CFF"))
ggplot(score.m[grepl("score", score.m$variable)&grepl("integral|peak_amp|speed",score.m$feature),], aes(variable, value))+
  geom_bar(stat= "identity", aes(fill = variable))+
  facet_wrap(~feature, nrow = 1)+ ylab("delta RSI")+
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.3,position=position_dodge(1))+
  # scale_fill_manual(values=colors_list)+
  theme_bw(base_size = 14)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#######################################################################
# Figure 7: Response Dynamics allows better polarization state identification----

#machine learning model LASSO on dynamics----
dynamics = rbind(cbind(dynamics_M0, type = "M0"),
                 cbind(dynamics_M1, type = "M1"),
                 cbind(dynamics_M2, type = "M2")) #M1, M2
dynamics = dynamics_M0M1M2 #use this, scaled across all polstates together instead of separately
head(dynamics)
collect_lasso = data.frame()
collect_f1 = data.frame()
collect_mi = read.delim("./infotheo/dynamics_k20_trajectoryfeatures_singlegene_SLEMI_polarization_specificity_perstim_scaledallM0M1M2.txt")
for (i in c("peak_amp","peak_amp_lfc", "integral","speed1hr")){
  for (s in c("CpG","IFNb", "LPS","P3CSK","PIC","TNF")){
    
    # i= "integral"
    # s= "LPS"
    print(i)
    print(s)
    
    collect_mi.subset = collect_mi[(collect_mi$feature==i)&(collect_mi$stimulus==s),]
    collect_mi.subset = collect_mi.subset[order(collect_mi.subset$cc, decreasing = T),]
    geneset = collect_mi.subset$gene[1:5] #get top 3 genes
    
    my.dataframe.subset = dynamics[dynamics$stimulus==s, c("cell","gene","stimulus","type", i)]
    my.dataframe.subset.dcast = dcast(my.dataframe.subset, cell+stimulus+type ~ gene, value.var = i)
    
    
    if(0){
      # ggplot(my.dataframe.subset.dcast, aes(type, Tgtp1))+geom_point(position = "jitter")
      # ggplot(my.dataframe.subset.dcast, aes(Fgl2, Tgtp1, color = type))+geom_point(position = "jitter")
      library(rgl);library(plot3Drgl);
      colors_list = c(M0= "dark red", M1= "#00BA38", M2= "#619CFF")
      with(my.dataframe.subset.dcast, scatter3D(x=Fgl2, y=Tgtp1, z=Irf1, colvar = as.integer(as.factor(type)),
                                                col = colors_list, bty = "b2",pch = 19, cex = 0.6, alpha = 0.3,   ticktype="detailed",
                                                xlab = "Fgl2", ylab = "Tgtp1", zlab = "Irf1"))
    }
    
    my.dataframe.subset.dcast = cbind(label = paste0(my.dataframe.subset.dcast$stimulus,"_",
                                                     my.dataframe.subset.dcast$type),my.dataframe.subset.dcast)
    
    library(caret);library(glmnet)
    set.seed(1)
    inTraining <- createDataPartition(my.dataframe.subset.dcast$label, p = .7, list = FALSE)
    training <- my.dataframe.subset.dcast[ inTraining,c("label", geneset)] #c("label", "Fgl2","Tgtp1","Irf1")
    testing  <- my.dataframe.subset.dcast[-inTraining,c("label", geneset)] #c("label", "Fgl2","Tgtp1","Irf1")
    
    #define response variable
    y <- as.factor(training$label)
    
    #define matrix of predictor variables  
    # x <- data.matrix(makeX(training[,-c(1:4)], na.impute = T))
    x <- data.matrix(makeX(training[,-1], na.impute = T)) #if selecting genes
    
    #perform k-fold cross-validation to find optimal lambda value
    cv_model <- cv.glmnet(x, y, alpha = 1, family="multinomial")
    plot(cv_model)
    
    #find optimal lambda value that minimizes test MSE
    best_lambda <- cv_model$lambda.min
    best_lambda
    
    #find coefficients of best model
    best_model <- glmnet(x, y, alpha = 1, family="multinomial",lambda = best_lambda)
    # plot(best_model, xvar = "dev", label = TRUE, type.coef = "coef")
    
    if(0){ #get betas and count
      tmp_coeffs = coef(best_model)
      beta <- Reduce(cbind, tmp_coeffs)
      beta <- beta[apply(beta != 0, 1, any),]
      colnames(beta) <- names(tmp_coeffs)
      beta 
      beta.frame = data.frame(beta)
      collect_lasso = rbind(collect_lasso, data.frame(feature = i, stimulus=s, num_vars=nrow(beta.frame)))
    }
    
    # split test-train
    if(1){ #get f1 scores on gene subsets
      #define new observation
      # new = data.matrix(testing[,-c(1:4)])
      new = data.matrix(testing[,-1])
      y_predicted = predict(best_model, s = best_lambda, newx = new)
      y.frame = data.frame(y_predicted[,,1])
      y.frame$predicted = colnames(y.frame)[apply(y.frame,1,which.max)]
      y.frame$actual = testing$label
      confusion = confusionMatrix(data = as.factor(y.frame$predicted), as.factor(y.frame$actual))
      values = as.data.frame(confusion$byClass)
      collect_f1 = rbind(collect_f1, data.frame(feature=i, stimulus = s, f1 = t(data.frame(values$F1))))
      
      # confusion.table = (confusion$table)
      # confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
      # p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
      #            # colorRampPalette(c("lightblue", "white", "red"))(50),
      #            # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
      #            breaks =  seq(0, 1, by = .01),
      #            color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
      #            display_numbers = T, fontsize_number = 14) 
      # 
      # y.frame$cell = seq(1:nrow(y.frame))
      # y.frame.m = melt(y.frame, id.vars = c("actual","cell"))
      # ggplot(y.frame.m, aes(variable, value, color = actual))+geom_violin()+
      #   geom_point(position = "jitter",size=0.1)+facet_grid(~actual)+
      #   # geom_line(aes(group = cell), size = 0.1)+
      #   theme_bw(base_size = 14)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    }
  }
}
# write.table(collect_lasso,"./infotheo/LASSO_dynamics_polarization_specificity_perstim.txt",sep="\t",quote=F,row.names = F)
# write.table(collect_lasso,"./infotheo/LASSO_dynamics_polarization_specificity_perstim_scaledallM0M1M2.txt",sep="\t",quote=F,row.names = F)
# write.table(collect_f1,"./infotheo/LASSO_dynamics_F1scoreTop5genes_polarization_specificity_perstim_scaledallM0M1M2.txt",sep="\t",quote=F,row.names = F)

#machine learning model LASSO on timepts----
collect = data.frame()
macro = readRDS(paste0("./output/macrophage_M0M1M2_combined_500genes_DBEC.rds"))
table(macro$stimulus)
collect_f1 = data.frame()
collect_mi = read.delim("./infotheo/timepoint0hr_singlegene_SLEMI_polarization_specificity_perstim.txt")

for (i in c("1hr","3hr","8hr")){
  # for (i in c("0.0hr")){
  print(i)
  macro.subset = subset(macro, subset = timept == i)
  table(macro.subset$stimulus)
  data = macro.subset[["ISnorm"]]@data
  data = data.frame(data)
  meta = macro.subset@meta.data
  colnames(data) = paste0(meta$stimulus, "_",meta$type)
  my.dataframe = cbind(label = colnames(data), data.frame(t(data)))
  head(my.dataframe)[1:5]
  table(my.dataframe$label)
  rownames(my.dataframe) = seq(1:nrow(my.dataframe))
  library(glmnet)
  
  for (s in c("CpG","IFNb", "LPS", "P3CSK","PIC","TNF")){
    # for (s in c("Unstim")){
    print(i)
    print(s)
    
    collect_mi.subset = collect_mi[(collect_mi$feature==i)&(collect_mi$stimulus==s),]
    collect_mi.subset = collect_mi.subset[order(collect_mi.subset$cc, decreasing = T),]
    geneset = collect_mi.subset$gene[1:5] #get top 3 genes
    
    my.dataframe.subset = my.dataframe[grepl(s, my.dataframe$label), ]
    
    if(0){
      ggplot(my.dataframe.subset, aes(label, Tgtp1))+geom_point(position = "jitter")
      ggplot(my.dataframe.subset, aes(Fgl2, Tgtp1, color = label))+geom_point(position = "jitter")
      library(rgl);library(plot3Drgl);
      colors_list = c(M0= "dark red", M1= "#00BA38", M2= "#619CFF")
      with(my.dataframe.subset, scatter3D(x=Fgl2, y=Tgtp1, z=Irf1, colvar = as.integer(as.factor(label)), 
                                          col = colors_list, bty = "b2",pch = 19, cex = 0.6, alpha = 0.3,   ticktype="detailed", 
                                          xlab = "Fgl2", ylab = "Tgtp1", zlab = "Irf1"))
    }
    
    library(caret)
    set.seed(1)
    inTraining <- createDataPartition(my.dataframe.subset$label, p = .7, list = FALSE)
    training <- my.dataframe.subset[ inTraining,c("label", geneset)] #Fgl2, Tgtp1, Irf1
    testing  <- my.dataframe.subset[-inTraining,c("label", geneset)]
    
    #define response variable
    y <- as.factor(training$label)
    
    #define matrix of predictor variables
    x <- data.matrix(training[,-1])
    
    #perform k-fold cross-validation to find optimal lambda value
    cv_model <- cv.glmnet(x, y, alpha = 1, family="multinomial")
    
    #find optimal lambda value that minimizes test MSE
    best_lambda <- cv_model$lambda.min
    best_lambda
    
    #find coefficients of best model
    best_model <- glmnet(x, y, alpha = 1, family="multinomial",lambda = best_lambda)
    # plot(best_model, xvar = "dev", label = TRUE, type.coef = "coef")
    
    if(0){ # for collect LASSO stats
      tmp_coeffs = coef(best_model)
      beta <- Reduce(cbind, tmp_coeffs)
      beta <- beta[apply(beta != 0, 1, any),]
      colnames(beta) <- names(tmp_coeffs)
      beta 
      beta.frame = data.frame(beta)
      collect = rbind(collect, data.frame(feature = i, stimulus=s, num_vars=nrow(beta.frame)))
    }
    ########################
    #use lasso regression model to predict response value
    # split test-train
    
    #define new observation
    new = data.matrix(testing[,-1])
    y_predicted = predict(best_model, s = best_lambda, newx = new)
    y.frame = data.frame(y_predicted[,,1])
    y.frame$predicted = colnames(y.frame)[apply(y.frame,1,which.max)]
    y.frame$actual = testing$label
    confusion = confusionMatrix(data = as.factor(y.frame$predicted), as.factor(y.frame$actual))
    confusion.table = (confusion$table)
    confusion.table = sweep(confusion.table, 2, colSums(confusion.table), FUN = '/')
    values = as.data.frame(confusion$byClass)
    collect_f1 = rbind(collect_f1, data.frame(feature=i, stimulus = s, f1 = t(data.frame(values$F1))))
    
    # p=pheatmap(confusion.table, cluster_rows = F, cluster_cols = F, 
    #            # colorRampPalette(c("lightblue", "white", "red"))(50),
    #            # gaps_col = c(5,10,15,20,25,30), gaps_row = c(5,10,15,20,25,30), 
    #            breaks =  seq(0, 1, by = .01),
    #            color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(seq(0, 1, by = .01) )),
    #            display_numbers = T, fontsize_number = 14) 
    
    # y.frame$cell = seq(1:nrow(y.frame))
    # y.frame.m = melt(y.frame, id.vars = c("actual","cell"))
    # ggplot(y.frame.m, aes(variable, value, color = actual))+geom_violin()+
    #   geom_point(position = "jitter",size=0.1)+facet_grid(~actual)+
    #   # geom_line(aes(group = cell), size = 0.1)+
    #   theme_bw(base_size = 14)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    
  }
} 
# write.table(collect, "./infotheo/LASSO_timepoints_polarization_specificity_perstim.txt",sep="\t",quote=F,row.names = F)
# write.table(collect_f1, "./infotheo/LASSO_timepoints_F1scoreTop5genes_polarization_specificity_perstim.txt",sep="\t",quote=F,row.names = F)

collect = read.delim("./infotheo/LASSO_timepoints_F1scoreTop5genes_polarization_specificity_perstim.txt")

#plot all LASSO----
collect_all = rbind(collect, collect_lasso)
colors_list = (c(Unstim="darkgray",CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3"))
ggplot(collect_all, aes(stimulus, num_vars)) +geom_bar(stat = "identity",aes(fill=stimulus))+
  facet_wrap(~feature, nrow = 2)+
  scale_fill_manual(values = colors_list)+ ylab("# LASSO-selected variables")+
  theme_bw(base_size = 14)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#plot all F1 scores----
collect = read.delim("./infotheo/LASSO_timepoints_F1scoreTop5genes_polarization_specificity_perstim.txt")
collect_f1=read.delim("./infotheo/LASSO_dynamics_F1scoreTop5genes_polarization_specificity_perstim_scaledallM0M1M2.txt")
collect_all = rbind(collect, collect_f1)
collect_all$f1.avg = rowMeans(collect_all[, c(3:5)])
colors_list = (c(Unstim="darkgray",CpG="#F8766D", IFNb="#B79F00",LPS= "#00BA38",P3CSK= "#00BFC4",PIC= "#619CFF",TNF= "#F564E3"))
ggplot(collect_all, aes(stimulus, f1.avg)) +geom_bar(stat = "identity",aes(fill=stimulus))+
  facet_wrap(~feature, nrow = 1)+
  scale_fill_manual(values = colors_list)+ ylab("F1 score - top 5 genes")+
  theme_bw(base_size = 14)+theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

