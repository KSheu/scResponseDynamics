#scResponseDynamics_script: imputation of single-cell response trajectories

source('trajectory4.R')

library(data.table);library(reshape2);library(pheatmap);
library(RColorBrewer);library(ggpubr);library(ggplot2);
library(ksheu.library1);library(SLEMI)library(Seurat);library(ggpubr);
setwd("F://scRNAseq_macro/scRNAseq_macro/")


###############################################################################
# Figure 1: Simulations of single-cell gene expression trajectories ----
###############################################################################
library(deSolve);library(optimx);library(gdata); library(reshape2) ; library(ggalt); library(dplyr);library(ggplot2)
library(dynamicTreeCut); library(gridExtra);library(ggpubr);library(amap); library(readxl);library(pheatmap)
library(RColorBrewer);library(reshape2);library(ggdendro);library(grid);library(ape);library(extrafont)
library(minpack.lm) # library for least squares fit using levenberg-marquart algorithm
library(rootSolve); library(FME) ; library(pso)#for fitting
library(greybox);library(truncnorm);library(ggalluvial);library(qlcMatrix)
library(torch);library(splines);library(smplot2)
setwd("F://scRNAseq_macro/scRNAseq_macro/")


# simulations of single cell gene trajectories prior to making data measurements

# bulk input trajectories (tables from Cheng et al. 2017) --------------------------------------------------
tf.ap1=read.delim("F://scRNAseq_macro/model_trajectory/2017_CellSystems_cheng_hoffmann/TableS2_AP1activationBMDM.txt")
tf.nfkb=read.delim("F://scRNAseq_macro/model_trajectory/2017_CellSystems_cheng_hoffmann/TableS2_NFkBactivationBMDM.txt")
tf.irf=read.delim("F://scRNAseq_macro/model_trajectory/2017_CellSystems_cheng_hoffmann/TableS2_IRFactivationBMDM.txt")
tf.p38=read.delim("F://scRNAseq_macro/model_trajectory/2017_CellSystems_cheng_hoffmann/Table_p38activationBMDM.txt")

tf.table.m = rbind(melt(data.frame(tf.ap1, tf = "AP1")), 
                   melt(data.frame(tf.nfkb, tf = "NFkB")), 
                   melt(data.frame(tf.irf, tf = "IRF")),
                   melt(data.frame(tf.p38, tf = "p38")))
tf.table.m$variable = as.numeric(gsub("X","",tf.table.m$variable))
colnames(tf.table.m) = c("stim", "tf","time", "value")
tf.table.m$stim = gsub("CPG", "CpG", tf.table.m$stim)
tf.table.m$stim = gsub("PAM", "P3CSK", tf.table.m$stim)

ggplot(tf.table.m, aes(time, value, group = tf))+geom_point(aes(color = tf), size = 2)+
  geom_line(aes(color = tf, group = tf), size = 1)+facet_wrap( ~stim, ncol = 7)+theme_bw(base_size = 16)
ggplot(tf.table.m, aes((time), value, group = tf))+geom_point(aes(color = tf), size = 2)+
  geom_line(aes(color = tf, group = tf), size = 0.5, linetype= "dashed")+
  geom_xspline(aes(color = tf, group = tf), spline_shape=-0.2, size=1)+
  facet_wrap( ~stim, ncol = 7)+theme_bw(base_size = 16) +xlim(0,300)

ggplot(tf.table.m[grepl("LPS|CpG|^PIC|P3CSK|TNF|IFNB", tf.table.m$stim),], aes((time), value, group = tf))+geom_point(aes(color = tf), size = 2)+
  geom_line(aes(color = tf, group = tf), size = 0.5, linetype= "dashed")+
  geom_xspline(aes(color = tf, group = tf), spline_shape=-0.2, size=1)+
  facet_wrap( ~stim, ncol = 3)+theme_bw(base_size = 16) +xlim(0,300)

# visualize interpolated inputs------------------------------------------------------
z = tf.table.m[tf.table.m$stim=="LPS"& tf.table.m$tf=="AP1", ]; plot(z$time, z$value)
input_lps.ap1 <-data.frame(xspline(z$time, z$value, shape=-0.3, draw=F)); ggplot(input_lps.ap1, aes(x,y)) +geom_point()+ylim(0,1)
z = tf.table.m[tf.table.m$stim=="LPS"& tf.table.m$tf=="NFkB", ]; plot(z$time, z$value)
input_lps.nfkb <-data.frame(xspline(z$time, z$value, shape=-0.3, draw=F)); ggplot(input_lps.nfkb, aes(x,y)) +geom_point()+ylim(0,1)
z = tf.table.m[tf.table.m$stim=="LPS"& tf.table.m$tf=="IRF", ]; plot(z$time, z$value)
input_lps.irf <-data.frame(xspline(z$time, z$value, shape= 0.4, draw=F)); ggplot(input_lps.irf, aes(x,y)) +geom_point()+ylim(0,1)

#generate single cell input functions based on bulk---------------------
collect_inputs_TFA_all = NULL
collect_inputs_TFN_all = NULL
collect_inputs_TFI_all = NULL
collect_inputs_p38_all = NULL

collect_Timepts_TFA_all = NULL
collect_Timepts_TFN_all = NULL
collect_Timepts_TFI_all = NULL
collect_Timepts_p38_all = NULL

cell_num = 100  
set.seed(123)
for (St_name in c("LPS","CpG","PIC","P3CSK","TNF","IFNB")){
  print(St_name)
  collect_inputs_TFA = data.frame(stimulus = c(rep(St_name,cell_num)))
  collect_inputs_TFN = data.frame(stimulus = c(rep(St_name,cell_num)))
  collect_inputs_TFI = data.frame(stimulus = c(rep(St_name,cell_num)))
  collect_inputs_p38 = data.frame(stimulus = c(rep(St_name,cell_num)))
  
  collect_Timepts_TFA = data.frame(stimulus = c(rep(St_name,cell_num)),Time0 = c(rep(0,cell_num)))
  collect_Timepts_TFN = data.frame(stimulus = c(rep(St_name,cell_num)),Time0 = c(rep(0,cell_num)))
  collect_Timepts_TFI = data.frame(stimulus = c(rep(St_name,cell_num)),Time0 = c(rep(0,cell_num)))
  collect_Timepts_p38 = data.frame(stimulus = c(rep(St_name,cell_num)),Time0 = c(rep(0,cell_num)))
  
  #AP1 profile
  TimepointsA = tf.table.m$time[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
  TFA_profile = tf.table.m$value[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
  
  #NFkB profile
  TimepointsN = tf.table.m$time[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
  TFN_profile = tf.table.m$value[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
  
  #IRF profile
  TimepointsI = tf.table.m$time[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
  TFI_profile = tf.table.m$value[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
  
  #p38 profile
  Timepointsp38 = tf.table.m$time[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
  p38_profile = tf.table.m$value[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
  
  for (datapt in seq(1:6)){
    print(datapt)
    
    # truncated normal distribution
    lnorm_ap1 <- rtruncnorm(cell_num, a=0, b=Inf,(TFA_profile[datapt]), log10(TFA_profile[datapt]+1))
    collect_inputs_TFA = cbind(collect_inputs_TFA,lnorm_ap1)
    
    lnorm_nfkb <- rtruncnorm(cell_num, a=0, b=Inf, (TFN_profile[datapt]), log10(TFN_profile[datapt]+1))
    collect_inputs_TFN = cbind(collect_inputs_TFN,lnorm_nfkb)
    
    if(datapt<6){
      lnorm_irf <- rtruncnorm(cell_num, a=0, b=Inf,(TFI_profile[datapt]), log10(TFI_profile[datapt]+1))
      collect_inputs_TFI = cbind(collect_inputs_TFI,lnorm_irf)
    }
    
    lnorm_p38 <- rtruncnorm(cell_num, a=0, b=Inf,(p38_profile[datapt]), log10(p38_profile[datapt]+1))
    collect_inputs_p38 = cbind(collect_inputs_p38,lnorm_p38)
    
  }
  for (datapt in c(2:6)){
    print(datapt)
    
    jitter_TimepointsA = rtruncnorm(cell_num,a=0, b=Inf,TimepointsA[datapt], log(TimepointsA[datapt]+1))
    collect_Timepts_TFA = cbind(collect_Timepts_TFA, jitter_TimepointsA)
    
    jitter_TimepointsN = rtruncnorm(cell_num,a=0, b=Inf,TimepointsN[datapt], log(TimepointsN[datapt]+1))
    collect_Timepts_TFN = cbind(collect_Timepts_TFN, jitter_TimepointsN)
    
    if(datapt<6){
      jitter_TimepointsI = rtruncnorm(cell_num,a=0, b=Inf,TimepointsI[datapt], log(TimepointsI[datapt]+1))
      collect_Timepts_TFI = cbind(collect_Timepts_TFI, jitter_TimepointsI)
    }
    
    jitter_Timepointsp38 = rtruncnorm(cell_num,a=0, b=Inf,Timepointsp38[datapt], log(Timepointsp38[datapt]+1))
    collect_Timepts_p38 = cbind(collect_Timepts_p38, jitter_Timepointsp38)
    
  }
  
  
  for (cell in seq(1:cell_num)){
    # print(cell)
    # TFA_profile_sc = collect_inputs_TFA[cell,-1]
    # TFN_profile_sc = collect_inputs_TFN[cell,-1]
    # TFI_profile_sc = collect_inputs_TFI[cell,-1]
    # p38_profile_sc = collect_inputs_p38[cell,-1]
    # 
    # ap1_sc <- approxfun(TimepointsA, TFA_profile_sc, rule =2)#(Time-tau)
    # nfkb_sc <- approxfun(TimepointsN, TFN_profile_sc, rule =2)#(Time-tau)
    # irf_sc <- approxfun(TimepointsI, TFI_profile_sc, rule =2)#(Time-tau)
    # p38_sc <- approxfun(Timepointsp38, p38_profile_sc, rule =2)#(Time)
    # 
    # curve(ap1_sc, -50,480, col = "darkorange",xlab = "time(mins)",  add = TRUE)
    # curve(nfkb_sc, -50, 480, col = "red", xlab = "time(mins)", add = TRUE)
    # curve(irf_sc, -50, 480, col = "darkgreen", xlab = "time(mins)",  add = TRUE)
    # curve(p38_sc, -50, 480, col = "blue", xlab = "time(mins)", add = TRUE)
  }
  
  collect_inputs_TFA_all = rbind(collect_inputs_TFA_all, collect_inputs_TFA)
  collect_inputs_TFN_all = rbind(collect_inputs_TFN_all, collect_inputs_TFN)
  collect_inputs_TFI_all = rbind(collect_inputs_TFI_all, collect_inputs_TFI)
  collect_inputs_p38_all = rbind(collect_inputs_p38_all, collect_inputs_p38)
  
  collect_Timepts_TFA_all = rbind(collect_Timepts_TFA_all, collect_Timepts_TFA)
  collect_Timepts_TFN_all = rbind(collect_Timepts_TFN_all, collect_Timepts_TFN)
  collect_Timepts_TFI_all = rbind(collect_Timepts_TFI_all, collect_Timepts_TFI)
  collect_Timepts_p38_all = rbind(collect_Timepts_p38_all, collect_Timepts_p38)
  
  
}

if (1){
  #replace negatives with 0 to avoid high baseline, not needed w/ truncated normal
  # collect_inputs_TFA_all[collect_inputs_TFA_all<0] <-0
  # collect_inputs_TFN_all[collect_inputs_TFN_all<0] <-0
  # collect_inputs_TFI_all[collect_inputs_TFI_all<0] <-0
  # collect_inputs_p38_all[collect_inputs_p38_all<0] <-0
  
  #normalize values to range 0-1
  collect_inputs_TFA_normalize = cbind(stimulus =collect_inputs_TFA_all$stimulus,
                                       (collect_inputs_TFA_all[,-1]- min(collect_inputs_TFA_all[,-1]))/(max(collect_inputs_TFA_all[,-1])-min(collect_inputs_TFA_all[,-1]))) #rescale each curveset 0-1 over all stims
  collect_inputs_TFN_normalize = cbind(stimulus =collect_inputs_TFN_all$stimulus,
                                       (collect_inputs_TFN_all[,-1]- min(collect_inputs_TFN_all[,-1]))/(max(collect_inputs_TFN_all[,-1])-min(collect_inputs_TFN_all[,-1]))) #rescale each curveset 0-1 over all stims
  collect_inputs_TFI_normalize = cbind(stimulus =collect_inputs_TFI_all$stimulus,
                                       (collect_inputs_TFI_all[,-1]- min(collect_inputs_TFI_all[,-1]))/(max(collect_inputs_TFI_all[,-1])-min(collect_inputs_TFI_all[,-1]))) #rescale each curveset 0-1 over all stims
  collect_inputs_p38_normalize = cbind(stimulus =collect_inputs_p38_all$stimulus,
                                       (collect_inputs_p38_all[,-1]- min(collect_inputs_p38_all[,-1]))/(max(collect_inputs_p38_all[,-1])-min(collect_inputs_p38_all[,-1]))) #rescale each curveset 0-1 over all stims
  
  
  par(mfrow = c(2, 3))
  for (St_name in c("CpG","IFNB","LPS","P3CSK","PIC","TNF")){
    print(St_name)
    #AP1 profile
    TimepointsA = tf.table.m$time[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
    TFA_profile = tf.table.m$value[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
    
    #NFkB profile
    TimepointsN = tf.table.m$time[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
    TFN_profile = tf.table.m$value[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
    
    #IRF profile
    TimepointsI = tf.table.m$time[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
    TFI_profile = tf.table.m$value[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
    
    #p38 profile
    Timepointsp38 = tf.table.m$time[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
    p38_profile = tf.table.m$value[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
    
    #plot TF activity forcing functions for the stimulus
    ap1 <- approxfun(TimepointsA, 
                     (TFA_profile-min(collect_inputs_TFA_all[,-1]))/(max(collect_inputs_TFA_all[,-1])-min(collect_inputs_TFA_all[,-1])), rule =2)#(Time-tau)
    nfkb <- approxfun(TimepointsN, 
                      (TFN_profile-min(collect_inputs_TFN_all[,-1]))/(max(collect_inputs_TFN_all[,-1])-min(collect_inputs_TFN_all[,-1])), rule =2)#(Time-tau)
    irf <- approxfun(TimepointsI, 
                     (TFI_profile-min(collect_inputs_TFI_all[,-1]))/(max(collect_inputs_TFI_all[,-1])-min(collect_inputs_TFI_all[,-1])), rule =2)#(Time-tau)
    p38 <- approxfun(Timepointsp38, 
                     (p38_profile-min(collect_inputs_p38_all[,-1]))/(max(collect_inputs_p38_all[,-1])-min(collect_inputs_p38_all[,-1])), rule =2)#(Time)
    curve(ap1, -10,480, col = "darkorange",xlab = "time(mins)", ylim =c(0,1.1),ylab="TF activity",main=St_name)
    curve(nfkb, -50, 480, col = "red", xlab = "time(mins)", add = TRUE)
    curve(irf, -50, 480, col = "darkgreen", xlab = "time(mins)",  add = TRUE)
    curve(p38, -50, 480, col = "blue", xlab = "time(mins)", add = TRUE)
    
    
    for (cell in seq(1:cell_num)){
      print(cell)
      TFA_profile_sc = collect_inputs_TFA_normalize[grepl(St_name, collect_inputs_TFA_normalize$stimulus),][cell,-1]
      TFN_profile_sc = collect_inputs_TFN_normalize[grepl(St_name, collect_inputs_TFN_normalize$stimulus),][cell,-1]
      TFI_profile_sc = collect_inputs_TFI_normalize[grepl(St_name, collect_inputs_TFI_normalize$stimulus),][cell,-1]
      p38_profile_sc = collect_inputs_p38_normalize[grepl(St_name, collect_inputs_p38_normalize$stimulus),][cell,-1]
      
      TimepointsA_sc = collect_Timepts_TFA_all[grepl(St_name, collect_Timepts_TFA_all$stimulus),][cell,-1]
      TimepointsN_sc = collect_Timepts_TFN_all[grepl(St_name, collect_Timepts_TFN_all$stimulus),][cell,-1]
      TimepointsI_sc = collect_Timepts_TFI_all[grepl(St_name, collect_Timepts_TFI_all$stimulus),][cell,-1]
      Timepointsp38_sc = collect_Timepts_p38_all[grepl(St_name, collect_Timepts_p38_all$stimulus),][cell,-1]
      
      # ap1_sc <- approxfun(TimepointsA_sc, TFA_profile_sc, rule =2)#(Time-tau)
      # nfkb_sc <- approxfun(TimepointsN_sc, TFN_profile_sc, rule =2)#(Time-tau)
      # irf_sc <- approxfun(TimepointsI_sc, TFI_profile_sc, rule =2)#(Time-tau)
      # p38_sc <- approxfun(Timepointsp38_sc, p38_profile_sc, rule =2)#(Time)
      # 
      # curve(ap1_sc, -50,480, col = alpha("darkorange",0.5),xlab = "time(mins)", add = TRUE)
      # curve(nfkb_sc, -50, 480, col = alpha("red",0.5), xlab = "time(mins)", add = TRUE)
      # curve(irf_sc, -50, 480, col = alpha("darkgreen",0.5), xlab = "time(mins)",  add = TRUE)
      # curve(p38_sc, -50, 480, col = alpha("blue",0.5), xlab = "time(mins)", add = TRUE)
      
      ap1_sc <- splinefun(TimepointsA_sc, TFA_profile_sc,method = "monoH.FC", ties = "constant")
      nfkb_sc <- splinefun(TimepointsN_sc, TFN_profile_sc,method = "monoH.FC", ties = "constant")
      irf_sc <- splinefun(TimepointsI_sc, TFI_profile_sc,method = "monoH.FC", ties = "constant")
      p38_sc <- splinefun(Timepointsp38_sc, p38_profile_sc,method = "monoH.FC", ties = "constant")
      
      curve(ap1_sc, -50,480, col = alpha("darkorange",0.1),xlab = "time(mins)", add = TRUE)
      curve(nfkb_sc, -50, 480, col = alpha("red",0.1), xlab = "time(mins)", add = TRUE)
      curve(irf_sc, -50, 480, col = alpha("darkgreen",0.1), xlab = "time(mins)",  add = TRUE)
      curve(p38_sc, -50, 480, col = alpha("blue",0.1), xlab = "time(mins)", add = TRUE)
      
    }
    
  }
  
  #plot with separated TFs
  par(mfrow = c(6, 4))
  for (St_name in c("CpG","IFNB","LPS","P3CSK","PIC","TNF")){
    print(St_name)
    #AP1 profile
    TimepointsA = tf.table.m$time[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
    TFA_profile = tf.table.m$value[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
    
    #NFkB profile
    TimepointsN = tf.table.m$time[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
    TFN_profile = tf.table.m$value[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
    
    #IRF profile
    TimepointsI = tf.table.m$time[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
    TFI_profile = tf.table.m$value[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
    
    #p38 profile
    Timepointsp38 = tf.table.m$time[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
    p38_profile = tf.table.m$value[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
    
    #plot TF activity forcing functions for the stimulus
    ap1 <- approxfun(TimepointsA, 
                     (TFA_profile-min(collect_inputs_TFA_all[,-1]))/(max(collect_inputs_TFA_all[,-1])-min(collect_inputs_TFA_all[,-1])), rule =2)#(Time-tau)
    nfkb <- approxfun(TimepointsN, 
                      (TFN_profile-min(collect_inputs_TFN_all[,-1]))/(max(collect_inputs_TFN_all[,-1])-min(collect_inputs_TFN_all[,-1])), rule =2)#(Time-tau)
    irf <- approxfun(TimepointsI, 
                     (TFI_profile-min(collect_inputs_TFI_all[,-1]))/(max(collect_inputs_TFI_all[,-1])-min(collect_inputs_TFI_all[,-1])), rule =2)#(Time-tau)
    p38 <- approxfun(Timepointsp38, 
                     (p38_profile-min(collect_inputs_p38_all[,-1]))/(max(collect_inputs_p38_all[,-1])-min(collect_inputs_p38_all[,-1])), rule =2)#(Time)
    
    curve(ap1, -50,480, col = "darkorange",xlab = "time(mins)", ylim =c(0,1),ylab="TF activity",main=paste0(St_name,"_AP1"))
    for (cell in seq(1:cell_num)){
      print(cell)
      TFA_profile_sc = collect_inputs_TFA_normalize[grepl(St_name, collect_inputs_TFA_normalize$stimulus),][cell,-1]
      TimepointsA_sc = collect_Timepts_TFA_all[grepl(St_name, collect_Timepts_TFA_all$stimulus),][cell,-1]
      # ap1_sc <- approxfun(TimepointsA_sc, TFA_profile_sc, rule =2)#(Time-tau)
      ap1_sc <- splinefun(TimepointsA_sc, TFA_profile_sc, method = "monoH.FC")#(Time-tau)
      curve(ap1_sc, -50,480, col = alpha("darkorange",0.05),xlab = "time(mins)", add = TRUE)
      
    }
    
    curve(nfkb, -50, 480, col = "red", xlab = "time(mins)",  ylim =c(0,1),ylab="TF activity",main=paste0(St_name,"_NFkB"))
    for (cell in seq(1:cell_num)){
      print(cell)
      TFN_profile_sc = collect_inputs_TFN_normalize[grepl(St_name, collect_inputs_TFN_normalize$stimulus),][cell,-1]
      TimepointsN_sc = collect_Timepts_TFN_all[grepl(St_name, collect_Timepts_TFN_all$stimulus),][cell,-1]
      # nfkb_sc <- approxfun(TimepointsN_sc, TFN_profile_sc, rule =2)#(Time-tau)
      nfkb_sc <- splinefun(TimepointsN_sc, TFN_profile_sc,  method = "monoH.FC")#(Time-tau)
      curve(nfkb_sc, -50, 480, col = alpha("red",0.05), xlab = "time(mins)", add = TRUE)
    }
    
    curve(irf, -50, 480, col = "darkgreen", xlab = "time(mins)",   ylim =c(0,1),ylab="TF activity",main=paste0(St_name,"_IRF"))
    for (cell in seq(1:cell_num)){
      print(cell)
      TFI_profile_sc = collect_inputs_TFI_normalize[grepl(St_name, collect_inputs_TFI_normalize$stimulus),][cell,-1]
      TimepointsI_sc = collect_Timepts_TFI_all[grepl(St_name, collect_Timepts_TFI_all$stimulus),][cell,-1]
      # irf_sc <- approxfun(TimepointsI_sc, TFI_profile_sc, rule =2)#(Time-tau)
      irf_sc <- splinefun(TimepointsI_sc, TFI_profile_sc, method = "monoH.FC" )#(Time-tau)
      curve(irf_sc, -50, 480, col = alpha("darkgreen",0.05), xlab = "time(mins)",  add = TRUE)
      
    }
    
    curve(p38, -50, 480, col = "blue", xlab = "time(mins)",  ylim =c(0,1),ylab="TF activity",main=paste0(St_name,"_p38"))
    for (cell in seq(1:cell_num)){
      print(cell)
      p38_profile_sc = collect_inputs_p38_normalize[grepl(St_name, collect_inputs_p38_normalize$stimulus),][cell,-1]
      Timepointsp38_sc = collect_Timepts_p38_all[grepl(St_name, collect_Timepts_p38_all$stimulus),][cell,-1]
      # p38_sc <- approxfun(Timepointsp38_sc, p38_profile_sc, rule =2)#(Time)
      p38_sc <- splinefun(Timepointsp38_sc, p38_profile_sc, method = "monoH.FC")#(Time)
      # p38_sc <- data.frame(xspline(Timepointsp38_sc, p38_profile_sc, shape=-0.3, draw=F))
      curve(p38_sc, -50, 480, col = alpha("blue",0.05), xlab = "time(mins)", add = TRUE)
      
    }
  }
  
}

#plot in inputs in heatmap form and write out dfs--------
if(1){ 
  
  plot_list=list()
  count=1
  for (St_name in c("CpG","IFNB","LPS","P3CSK","PIC","TNF")){
    
    print(St_name)
    heatmapdf_ap1 = NULL
    heatmapdf_nfkb = NULL
    heatmapdf_irf = NULL
    heatmapdf_p38 = NULL
    for (cell in seq(1:cell_num)){
      print(cell)
      TFA_profile_sc = collect_inputs_TFA_normalize[grepl(St_name, collect_inputs_TFA_normalize$stimulus),][cell,-1]
      TFN_profile_sc = collect_inputs_TFN_normalize[grepl(St_name, collect_inputs_TFN_normalize$stimulus),][cell,-1]
      TFI_profile_sc = collect_inputs_TFI_normalize[grepl(St_name, collect_inputs_TFI_normalize$stimulus),][cell,-1]
      p38_profile_sc = collect_inputs_p38_normalize[grepl(St_name, collect_inputs_p38_normalize$stimulus),][cell,-1]
      
      TimepointsA_sc = collect_Timepts_TFA_all[grepl(St_name, collect_Timepts_TFA_all$stimulus),][cell,-1]
      TimepointsN_sc = collect_Timepts_TFN_all[grepl(St_name, collect_Timepts_TFN_all$stimulus),][cell,-1]
      TimepointsI_sc = collect_Timepts_TFI_all[grepl(St_name, collect_Timepts_TFI_all$stimulus),][cell,-1]
      Timepointsp38_sc = collect_Timepts_p38_all[grepl(St_name, collect_Timepts_p38_all$stimulus),][cell,-1]
      
      # ap1_sc <- approxfun(TimepointsA_sc, TFA_profile_sc, rule =2)#(Time-tau)
      # nfkb_sc <- approxfun(TimepointsN_sc, TFN_profile_sc, rule =2)#(Time-tau)
      # irf_sc <- approxfun(TimepointsI_sc, TFI_profile_sc, rule =2)#(Time-tau)
      # p38_sc <- approxfun(Timepointsp38_sc, p38_profile_sc, rule =2)#(Time)
      
      ap1_sc <- splinefun(TimepointsA_sc, TFA_profile_sc,method = "monoH.FC")
      nfkb_sc <- splinefun(TimepointsN_sc, TFN_profile_sc,method = "monoH.FC")
      irf_sc <- splinefun(TimepointsI_sc, TFI_profile_sc,method = "monoH.FC")
      p38_sc <- splinefun(Timepointsp38_sc, p38_profile_sc,method = "monoH.FC")
      
      x = seq(-50, 480, length.out = 530)
      curve_ap1 = ap1_sc(x)
      curve_nfkb= nfkb_sc(x)
      curve_irf= irf_sc(x)
      curve_p38=p38_sc(x)
      
      heatmapdf_ap1 = rbind(heatmapdf_ap1, curve_ap1)
      heatmapdf_nfkb = rbind(heatmapdf_nfkb, curve_nfkb)
      heatmapdf_irf = rbind(heatmapdf_irf, curve_irf)
      heatmapdf_p38 = rbind(heatmapdf_p38, curve_p38)
      
    }
    
    heatmapdf_ap1= data.frame(heatmapdf_ap1)
    heatmapdf_nfkb= data.frame(heatmapdf_nfkb)
    heatmapdf_irf= data.frame(heatmapdf_irf)
    heatmapdf_p38= data.frame(heatmapdf_p38)
    
    heatmapdf_ap1[heatmapdf_ap1<0] <-0
    heatmapdf_nfkb[heatmapdf_nfkb<0] <-0
    heatmapdf_irf[heatmapdf_irf<0] <-0
    heatmapdf_p38[heatmapdf_p38<0] <-0
    
    
    #write out the dataframes
    # write.table(heatmapdf_ap1,paste0("./generated_scTFinputs/scTFinputs_",St_name,"_AP1.txt"),quote = F,sep="\t",row.names = F)
    # write.table(heatmapdf_nfkb,paste0("./generated_scTFinputs/scTFinputs_",St_name,"_NFkB.txt"),quote = F,sep="\t",row.names = F)
    # write.table(heatmapdf_irf,paste0("./generated_scTFinputs/scTFinputs_",St_name,"_IRF.txt"),quote = F,sep="\t",row.names = F)
    # write.table(heatmapdf_p38,paste0("./generated_scTFinputs/scTFinputs_",St_name,"_p38.txt"),quote = F,sep="\t",row.names = F)
    
    p1=pheatmap(heatmapdf_ap1[order(apply(heatmapdf_ap1[,50:170],1,max), decreasing = T),], scale="none",
                cluster_cols = F, cluster_rows = F,border_color = NA,
                color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                breaks = c(0,seq(0.01,0.5,length=100),1),
                show_colnames = F, show_rownames = F, main = paste0(St_name,"_AP1"))
    plot_list[[count]] =p1[[4]]; count=count+1
    p2=pheatmap(heatmapdf_nfkb[order(apply(heatmapdf_nfkb[,50:170],1,max), decreasing = T),], scale="none",
                cluster_cols = F, cluster_rows = F,border_color = NA,
                color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                breaks = c(0,seq(0.01,0.5,length=100),1),
                show_colnames = F, show_rownames = F, main = paste0(St_name,"_NFkB"))
    plot_list[[count]] =p2[[4]]; count=count+1
    p3=pheatmap(heatmapdf_irf[order(apply(heatmapdf_irf[,50:170],1,max), decreasing = T),], scale="none",
                cluster_cols = F, cluster_rows = F,border_color = NA,
                color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                breaks = c(0,seq(0.01,0.5,length=100),1),
                show_colnames = F, show_rownames = F, main = paste0(St_name,"_IRF"))
    plot_list[[count]] =p3[[4]]; count=count+1
    p4=pheatmap(heatmapdf_p38[order(apply(heatmapdf_p38[,50:170],1,max), decreasing = T),], scale="none",
                cluster_cols = F, cluster_rows = F,border_color = NA,
                color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                breaks = c(0,seq(0.01,0.5,length=100),1),
                show_colnames = F, show_rownames = F, main = paste0(St_name,"_p38"))
    plot_list[[count]] =p4[[4]]; count=count+1
    
    
  }
  grid.arrange(grobs=plot_list, nrow=4)
}
#check p38 curve effect on kdeg----
p38=read.delim("./generated_scTFinputs/scTFinputs_LPS_p38.txt")
kdeg = matrix(log(2)/30, nrow = nrow(p38),ncol = ncol(p38))
kdeg_p38=log(2)/((log(2)/kdeg) + 480*p38)
pheatmap(kdeg_p38, cluster_cols = F, cluster_rows = T, show_colnames = F, main="kdeg_p38")


# differential equations----------------------------------------------
TFA_profile_df_all = NULL
TFN_profile_df_all = NULL
TFI_profile_df_all = NULL
p38_profile_df_all = NULL

for (s in c("CpG", "IFNB", "LPS","P3CSK", "PIC", "TNF")){
  TFA_tmp = read.delim(paste0("./generated_scTFinputs/scTFinputs_",s,"_AP1.txt"))
  TFN_tmp = read.delim(paste0("./generated_scTFinputs/scTFinputs_",s,"_NFkB.txt"))
  TFI_tmp = read.delim(paste0("./generated_scTFinputs/scTFinputs_",s,"_IRF.txt"))
  p38_tmp = read.delim(paste0("./generated_scTFinputs/scTFinputs_",s,"_p38.txt"))

  TFA_profile_df_all = rbind(TFA_profile_df_all, data.frame(stim=s,TFA_tmp))
  TFN_profile_df_all = rbind(TFN_profile_df_all, data.frame(stim=s,TFN_tmp))
  TFI_profile_df_all = rbind(TFI_profile_df_all, data.frame(stim=s,TFI_tmp))
  p38_profile_df_all = rbind(p38_profile_df_all, data.frame(stim=s,p38_tmp))
  
}

odeModel_steadystate <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)),{
    
    St_name<-stims[Pars[[c("stimulus")]] ]
    print(paste0("currently simulating Stim: ",St_name))
    # St_name = "TNF"
    
    gene.clust = GRS_list[Pars[c("GRS")]  ]
    # gene.clust = collect_cost_all$best[collect_cost_all$gene==gene]
    # gene.clust = "G2L"#genetypes[Pars[[c("genetypes")]] ]  #mC testing just one GRS for now
    
    # #AP1 profile
    # TimepointsA = tf.table.m$time[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
    # TFA_profile = tf.table.m$value[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
    # 
    # #NFkB profile
    # TimepointsN = tf.table.m$time[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
    # TFN_profile = tf.table.m$value[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
    # 
    # #IRF profile
    # TimepointsI = tf.table.m$time[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
    # TFI_profile = tf.table.m$value[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
    # 
    # #p38 profile
    # Timepointsp38 = tf.table.m$time[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
    # p38_profile = tf.table.m$value[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
    
 
    n = Pars[c("n")]
    ktA<-Pars[c("ktA")]
    ktN<-Pars[c("ktN")]
    ktI<-Pars[c("ktI")]
    k0 <- Pars[c("k0")]
    k_deg = Pars[c("k_deg")]
    k_syn = k_deg 
    # tau = Pars[c("tau")]
    tau = 0
    cell = Pars[c("cell")]
    
    #TF activity forcing function for one input
    # ap1 <- approxfun(TimepointsA, TFA_profile, rule =2)(Time-tau)
    # nfkb <- approxfun(TimepointsN, TFN_profile, rule =2)(Time-tau)
    # irf <- approxfun(TimepointsI, TFI_profile, rule =2)(Time-tau)
    # p38 <- approxfun(Timepointsp38, p38_profile, rule =2)(Time)
    # curve(ap1, -50,480, col = "darkorange",xlab = "time(mins)",ylim = c(0,1))
    # curve(nfkb, -50, 480, col = "red", xlab = "time(mins)", add = TRUE)
    # curve(irf, -50, 480, col = "darkgreen", xlab = "time(mins)",  add = TRUE)
    # curve(p38, -50, 480, col = "blue", xlab = "time(mins)", add = TRUE)
    
    
    # TFA_profile_df = read.delim(paste0("./generated_scTFinputs/scTFinputs_",St_name,"_AP1.txt"))
    # TFN_profile_df = read.delim(paste0("./generated_scTFinputs/scTFinputs_",St_name,"_NFkB.txt"))
    # TFI_profile_df = read.delim(paste0("./generated_scTFinputs/scTFinputs_",St_name,"_IRF.txt"))
    # p38_profile_df = read.delim(paste0("./generated_scTFinputs/scTFinputs_",St_name,"_p38.txt"))
    
    TFA_profile_df = TFA_profile_df_all[TFA_profile_df_all$stim==St_name,-1]
    TFN_profile_df = TFN_profile_df_all[TFN_profile_df_all$stim==St_name,-1]
    TFI_profile_df = TFI_profile_df_all[TFI_profile_df_all$stim==St_name,-1]
    p38_profile_df = p38_profile_df_all[p38_profile_df_all$stim==St_name,-1]
    
    TimeRange = c(-49:480)
    
    #TF activity forcing function for scTFinputs
    ap1 <- approxfun(TimeRange, TFA_profile_df[cell,], rule =2)(Time-tau)
    nfkb <- approxfun(TimeRange, TFN_profile_df[cell,], rule =2)(Time-tau)
    irf <- approxfun(TimeRange, TFI_profile_df[cell,], rule =2)(Time-tau)
    p38 <- approxfun(TimeRange, p38_profile_df[cell,], rule =2)(Time)
    
    
    
    # #single gates
    # fA<-(1.0-k0)*((ktA*ap1)^n/(1.0+ktA*ap1)^n)+k0
    # fN<-(1.0-k0)*((ktN*nfkb)^n/(1.0+ktN*nfkb)^n)+k0
    # fI<-(1.0-k0)*((ktI*irf)^n/(1.0+ktI*irf)^n)+k0
    # 
    # #OR gate
    # fIN<-(1.0-k0)*((ktN*nfkb+ktI*irf+ktN*ktI*nfkb*irf)^n/(1.0+ktN*nfkb+ktI*irf+ktN*ktI*nfkb*irf)^n)+k0
    
    
    
    #single gates
    fA<-(1.0-k0)*((ktA*ap1)^n/(1.0+(ktA*ap1)^n))+k0
    fN<-(1.0-k0)*((ktN*nfkb)^n/(1.0+(ktN*nfkb)^n))+k0
    fI<-(1.0-k0)*((ktI*irf)^n/(1.0+(ktI*irf)^n))+k0
    
    #OR gate
    fIN<-(1.0-k0)*(((ktN*nfkb)^n+(ktI*irf)^n+(ktN*ktI*nfkb*irf)^n)/(1.0+(ktN*nfkb)^n+(ktI*irf)^n+(ktN*ktI*nfkb*irf)^n))+k0
    
    #AND gate
    fAIN<-(1.0-k0)*((ktA*ktN*ktI*ap1*nfkb*irf)^n /(1.0+ktA*ap1+ktN*nfkb+ktI*irf+ktA*ktN*ap1*nfkb+
                                                     ktN*ktI*nfkb*irf+ktA*ktI*ap1*irf+
                                                     ktA*ktN*ktI*ap1*nfkb*irf)^n)+k0
    # mRNA ODEs
    
    if (gene.clust =="AP1"){ #keeping all original names, but using modv3,p38input
      dmRNA <-k_syn*fA - k_deg*mRNA #mA, G1S
    }
    if (gene.clust =="NFkB"){
      dmRNA <-k_syn*fN - k_deg*mRNA #mC, G2L, #mB, G2S
    }
    if (gene.clust =="NFkB|p38"){
      dmRNA <-k_syn*fN - k_deg*mRNA  #mD, G10L
    }
    if (gene.clust =="NFkB|IRF"){
      dmRNA <-k_syn*fIN - k_deg*mRNA #mE, G7S
    }
    if (gene.clust =="IRF"){
      dmRNA <-k_syn*fI - k_deg*mRNA #mF, G3S,#mG, G3L
    }
    
    return(list(c(dmRNA) ))
  })
}
odeModel <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)),{
    
    St_name<-stims[Pars[[c("stimulus")]] ]
    # print(paste0("current St_name ",St_name))
    # St_name = "TNF"
    
    gene.clust = GRS_list[Pars[c("GRS")] ] 
    # gene.clust = collect_cost_all$best[collect_cost_all$gene==gene]
    # gene.clust = "G2L" #mC testing just one GRS for now
    
    # #AP1 profile
    # TimepointsA = tf.table.m$time[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
    # TFA_profile = tf.table.m$value[tf.table.m$tf =="AP1" & tf.table.m$stim==St_name]
    # 
    # #NFkB profile
    # TimepointsN = tf.table.m$time[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
    # TFN_profile = tf.table.m$value[tf.table.m$tf =="NFkB" & tf.table.m$stim==St_name]
    # 
    # #IRF profile
    # TimepointsI = tf.table.m$time[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
    # TFI_profile = tf.table.m$value[tf.table.m$tf =="IRF" & tf.table.m$stim==St_name]
    # 
    # #p38 profile
    # Timepointsp38 = tf.table.m$time[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
    # p38_profile = tf.table.m$value[tf.table.m$tf =="p38" & tf.table.m$stim==St_name]
    
    n = Pars[c("n")]
    ktA<-Pars[c("ktA")]
    ktN<-Pars[c("ktN")]
    ktI<-Pars[c("ktI")]
    k0 <-Pars[c("k0")]
    k_deg = Pars[c("k_deg")]
    k_deg_p38 = Pars[c("k_deg")]
    # k_p38 = 1 #modifier of p38 effect on t1/2 #Pars[c("k_p38")]
    k_syn = k_deg 
    tau = Pars[c("tau")]
    cell = Pars[c("cell")]
    
    #TF activity forcing function 
    # ap1 <- approxfun(TimepointsA, TFA_profile, rule =2)(Time-tau)
    # nfkb <- approxfun(TimepointsN, TFN_profile, rule =2)(Time-tau)
    # irf <- approxfun(TimepointsI, TFI_profile, rule =2)(Time-tau)
    # p38 <- approxfun(Timepointsp38, p38_profile, rule =2)(Time)
    
    # TFA_profile_df = read.delim(paste0("./generated_scTFinputs/scTFinputs_",St_name,"_AP1.txt"))
    # TFN_profile_df = read.delim(paste0("./generated_scTFinputs/scTFinputs_",St_name,"_NFkB.txt"))
    # TFI_profile_df = read.delim(paste0("./generated_scTFinputs/scTFinputs_",St_name,"_IRF.txt"))
    # p38_profile_df = read.delim(paste0("./generated_scTFinputs/scTFinputs_",St_name,"_p38.txt"))
    
    TFA_profile_df = TFA_profile_df_all[TFA_profile_df_all$stim==St_name,-1]
    TFN_profile_df = TFN_profile_df_all[TFN_profile_df_all$stim==St_name,-1]
    TFI_profile_df = TFI_profile_df_all[TFI_profile_df_all$stim==St_name,-1]
    p38_profile_df = p38_profile_df_all[p38_profile_df_all$stim==St_name,-1]
    
    TimeRange = c(-49:480)
    
    #TF activity forcing function for scTFinputs
    ap1 <- approxfun(TimeRange, TFA_profile_df[cell,], rule =2)(Time-tau)
    nfkb <- approxfun(TimeRange, TFN_profile_df[cell,], rule =2)(Time-tau)
    irf <- approxfun(TimeRange, TFI_profile_df[cell,], rule =2)(Time-tau)
    p38 <- approxfun(TimeRange, p38_profile_df[cell,], rule =2)(Time)
    
    
    
    if (St_name=="LPS"|St_name=="P3CSK"|St_name=="CpG"){
      mRNA_life_mod = (log(2)/k_deg)+(480*p38)
      # print(mRNA_life_mod)
      k_deg_p38<-log(2)/mRNA_life_mod
    }
    
    if(0){
      #plot p38 vs kdeg
      curve(log(2)/((30)+(480*x)), 0,1, ylab="kdeg",xlab = "[p38]")
      curve(log(2)/((120)+(480*x)), 0,1, add = T, col = 'purple')
      curve(log(2)/((180)+(480*x)), 0,1, add = T, col = "gray")
      curve(log(2)/((300)+(480*x)), 0,1, add = T, col = "pink")
      
      #plot TFs vs ksyn
      curve((1.0-0.005)*((5.2*x)^3/(1.0+(5.2*x)^3))+0.005, 0,1, col = 'blue',ylab="synthesis",xlab = "[TF]")
      curve((1.0-0.005)*((0.32*x)^3/(1.0+(0.32*x)^3))+0.005, 0,1, add=T,col = 'blue',ylab="synthesis",xlab = "[TF]")
      curve((1.0-0.005)*((0.65*x)^3/(1.0+(0.65*x)^3))+0.005, 0,1, add = T, col = 'blue',ylab="synthesis",xlab = "[TF]")
      curve((1.0-0.005)*((1.3*x)^3/(1.0+(1.3*x)^3))+0.005, 0,1, add = T, col = 'blue',ylab="synthesis",xlab = "[TF]")
      curve((1.0-0.005)*((2.6*x)^3/(1.0+(2.6*x)^3))+0.005, 0,1, add = T,col = 'blue',ylab="synthesis",xlab = "[TF]")
      
      curve((1.0-0.005)*((5.2*x)^1/(1.0+(5.2*x)^1))+0.005, 0,1, add = T,col = 'skyblue',ylab="synthesis",xlab = "[TF]")
      curve((1.0-0.005)*((0.32*x)^1/(1.0+(0.32*x)^1))+0.005, 0,1, add=T,col = 'skyblue',ylab="synthesis",xlab = "[TF]")
      curve((1.0-0.005)*((0.65*x)^1/(1.0+(0.65*x)^1))+0.005, 0,1, add = T, col = 'skyblue',ylab="synthesis",xlab = "[TF]")
      curve((1.0-0.005)*((1.3*x)^1/(1.0+(1.3*x)^1))+0.005, 0,1, add = T, col = 'skyblue',ylab="synthesis",xlab = "[TF]")
      curve((1.0-0.005)*((2.6*x)^1/(1.0+(2.6*x)^1))+0.005, 0,1, add = T,col = 'skyblue',ylab="synthesis",xlab = "[TF]")
    }
    
    #single gates
    fA<-(1.0-k0)*((ktA*ap1)^n/(1.0+(ktA*ap1)^n))+k0
    fN<-(1.0-k0)*((ktN*nfkb)^n/(1.0+(ktN*nfkb)^n))+k0
    fI<-(1.0-k0)*((ktI*irf)^n/(1.0+(ktI*irf)^n))+k0
    
    #OR gate
    fIN<-(1.0-k0)*(((ktN*nfkb)^n+(ktI*irf)^n+(ktN*ktI*nfkb*irf)^n)/(1.0+(ktN*nfkb)^n+(ktI*irf)^n+(ktN*ktI*nfkb*irf)^n))+k0
    
    #AND gate
    fAIN<-(1.0-k0)*((ktA*ktN*ktI*ap1*nfkb*irf)^n /(1.0+ktA*ap1+ktN*nfkb+ktI*irf+ktA*ktN*ap1*nfkb+
                                                     ktN*ktI*nfkb*irf+ktA*ktI*ap1*irf+
                                                     ktA*ktN*ktI*ap1*nfkb*irf)^n)+k0
    # mRNA ODEs
    
    if (gene.clust =="AP1"){ #keeping all original names, but using modv3,p38input
      dmRNA <-k_syn*fA - k_deg*mRNA #mA, G1S
    }
    if (gene.clust =="NFkB"){
      dmRNA <-k_syn*fN - k_deg*mRNA #mC, G2L, #mB, G2S
    }
    if (gene.clust =="NFkB|p38"){
      dmRNA <-k_syn*fN - k_deg_p38*mRNA  #mD, G10L
    }
    if (gene.clust =="NFkB|IRF"){
      dmRNA <-k_syn*fIN - k_deg*mRNA #mE, G7S
    }
    if (gene.clust =="IRF"){
      dmRNA <-k_syn*fI - k_deg*mRNA #mF, G3S,#mG, G3L
    }
    
    return(list(c(dmRNA) ))
  })
}




# params and initialization-------------------------------------------------------------
pars = c(mRNA = 0.01, n = 1, ktA =0.48 , ktN = 1.3, ktI = 1.25, k0=.005,
         k_deg=log(2)/30, tau = 0.001, stimulus = NA, cell=NA, GRS = NA)
# pars = c(mRNA = 0.01, n = 1, ktA =0.5 , ktN = 0.5, ktI = 0.5, k0=.005,  
#          k_deg=log(2)/90, tau = 0.001, stimulus = NA, cell=NA, GRS = NA)

stims = "LPS" #c("CpG", "IFNB", "LPS","P3CSK", "PIC", "TNF")
stims = c("CpG", "IFNB", "LPS","P3CSK", "PIC", "TNF")
GRS_list = c("AP1","NFkB","NFkB|p38","NFkB|IRF","IRF")
cell_num=50
set.seed(123)
cells = sample(c(1:100),cell_num,replace = F)

solve_model_steadystate <- function(pars, times = seq(-300, 0, length=10)) {
  # times = seq(0, 480, length=500)
  cnt = 1
  for (i in stims){
    
    pars[c("stimulus")] = c(stimulus = which(stims==i) )
    # print(pars[c("stimulus")])
    
    for(j in cells){
      pars[c("cell")] = j
      
      for (k in GRS_list){
        pars[c("GRS")] = c(GRS = which(GRS_list==k) )
        
        ss.init <- pars[c("mRNA")]
        params <- pars[c("n", "ktA", "ktN", "ktI", "k0", "k_deg", "tau", "stimulus", "cell", "GRS")]
        out <- ode(ss.init, times, odeModel_steadystate, params, method = "lsode", atol = 1e-3, rtol = 1e-3)
        out.frame = as.data.frame(out)
        out.frame$label = i
        out.frame$cell = paste0("cell_",j )
        out.frame$GRS = k
        out.frame$params = list(pars[c("n", "ktA", "ktN", "ktI", "k0", "k_deg", "tau")])
        
        if(cnt==1){
          collect_out <- out.frame
          
        }else{
          collect_out <-rbind(collect_out, out.frame)
        }
        cnt=cnt+1
      }
    }
  }
  return(collect_out)
}


solve_model <- function(pars, times = seq(0, 480, length=300)) {
  # times = seq(0, 480, length=500)
  
  solve.ss = solve_model_steadystate(pars)
  steady.state = (solve.ss[nrow(solve.ss),2])
  
  cnt = 1
  for (i in stims){
    
    pars[c("mRNA")] <- steady.state
    
    # pars[c("mRNA")] <- runsteady(y = c(mRNA=0), #init
    #                              fun = odeModel_steadystate,
    #                              parms = pars, time = c(0,1e5))$y
    
    
    pars[c("stimulus")] = c(stimulus = which(stims==i) )
    print(i)
    
    for(j in cells){
      pars[c("cell")] = j
      print(paste0("cell_",j ))
      
      for (k in GRS_list){
        print(k)
        pars[c("GRS")] = c(GRS = which(GRS_list==k) )
        
        ss.init <- pars[c("mRNA")]
        params <- pars[c("n", "ktA", "ktN", "ktI", "k0",  "k_deg", "tau", "stimulus", "cell", "GRS")]
        out <- ode(ss.init, times, odeModel, params, method = "lsode",atol = 1e-3, rtol = 1e-3)
        out.frame = as.data.frame(out)
        out.frame$label = i
        out.frame$cell = paste0("cell_",j )
        out.frame$GRS = k
        out.frame$params = list(pars[c("n", "ktA", "ktN", "ktI", "k0", "k_deg", "tau")])
        
        if(cnt==1){
          collect_out <- out.frame
          
        }else{
          collect_out <-rbind(collect_out, out.frame)
        }
        cnt=cnt+1
      }
    }
  }
  
  return(collect_out)
}

# simulate model with one set of gene reg params-------------------------------------
par(mfrow=c(1,1))
solve = solve_model(pars)

solve = readRDS(paste0("./generated_scTFinputs/output_test_5GRS.rds"))

solve.truncate = solve[solve$time >-100, ]
ggplot(data = solve.truncate, aes(time, mRNA, color = cell))+
  geom_line(aes(group=cell), size = 0.5, alpha=0.9)+
  xlim(-10, 500)+ theme_bw(base_size = 10)+
  # stat_summary(aes(time, mRNA),geom="line", color = "black",fun = "mean", size=.3, linetype="dashed")+
  facet_grid(GRS~label, scales = "free")


# generate table of param sets to vary gene reg param--------
pars = c(mRNA = 0.01, n = 1, ktA =0.48 , ktN = 1.3, ktI = 1.25, k0=.005,
         k_deg=log(2)/30, tau = 0.001, stimulus = NA, cell=NA, GRS = NA)


param_sets = NULL
for (i in c(1,3)){
  for (j in c(0.25, 0.5, 1, 2,4)){
    for (k in c(30,120,180,300)){
      pars_new = c(mRNA = 0.01, n = i, 
                   ktA =j*pars[[c('ktA')]] , ktN = j*pars[[c('ktN')]], ktI = j*pars[[c('ktI')]], k0=.005,  
                   k_deg=(log(2)/k), tau = 0.001, stimulus = NA, cell=NA, GRS = NA)
      param_sets = rbind(param_sets, pars_new)
    }
    
  }
  
}
param_sets = data.frame(param_sets)

library(doParallel)
cl <- makeCluster(detectCores()-2)
registerDoParallel(cl)
#.export = c("odeModel", "odeModel_steadystate","solve_model","solve_model_steadystate")
foreach::foreach(i=seq(1:nrow(param_sets)), .packages=c("deSolve")) %dopar% {
  # .GlobalEnv$mRNA <- mRNA
  print(i)
  pars = unlist(param_sets[i,])
  solve = solve_model(pars)
  # saveRDS(solve,file = paste0("./generated_scTFinputs/output_test_5GRS_paramset",i,"_LPS.rds"))
}
stopCluster(cl)


# read back in the outputs and plot-----
solve.truncate = NULL
for (i in seq(1:nrow(param_sets))){
  print(i) #i=9 is average fit
  solve = readRDS(paste0("./generated_scTFinputs/output_test_5GRS_paramset",i,"_LPS.rds"))
  
  #plot for Fig1 (i=13)
  solve$GRS = factor(solve$GRS, levels = c("AP1","NFkB","NFkB|p38","NFkB|IRF","IRF"))
  ggplot(data = solve, aes(time/60, mRNA, color = GRS))+ xlab("time (hrs)")+xlim(0,8)+ 
    scale_color_manual(values = c(AP1="darkorange", IRF="forestgreen",NFkB= "red", `NFkB|p38`= "blue1",`NFkB|IRF`= "darkred"))+
    geom_line(aes(group=cell), size = 0.1, alpha=0.8)+
    theme_bw(base_size = 10)+ theme(legend.position = "None")+
    stat_summary(aes(time/60, mRNA),geom="line", color = "black",fun = "mean", linewidth=.3, alpha=0.4,linetype="dashed")+
    facet_wrap(~GRS, nrow=1,scales = "free")+ggtitle(solve$params[i])
  #plot for Fig1, cell_23, cell_50, cell_7, cell_32
  ggplot(data = solve[grepl("cell_32$",solve$cell)&grepl("NFkB$|AP1$|^IRF$",solve$GRS),], aes(time, mRNA, color = GRS))+
    scale_color_manual(values = c(AP1="darkorange", IRF="forestgreen",NFkB= "red", `NFkB|p38`= "blue1",`NFkB|IRF`= "darkred"))+
    geom_line(aes(group=GRS), size = 1, alpha=0.9)+
    xlim(-10, 500)+ theme_bw(base_size = 10)+theme_classic(base_size = 14)+theme(legend.position = "None")+
    ggtitle(solve$params[i])
  
  
  ggplot(data = solve, aes(time, mRNA, color = cell))+
    geom_line(aes(group=cell), size = 0.1, alpha=0.9)+
    xlim(-10, 500)+ theme_bw(base_size = 10)+ theme(legend.position = "None")+
    stat_summary(aes(time, mRNA),geom="line", color = "black",fun = "mean", linewidth=.3, linetype="dashed")+
    facet_grid(GRS~label, scales = "free")+ggtitle(solve$params[i])
  # print(p)
  ggsave(paste0("./generated_scTFinputs/output_test_5GRS_paramset",i,"_LPS.pdf"),
         width = 3.5,height = 7,units = "in")
  # solve.truncate = rbind(solve.truncate, solve)
}


#################################----
# rearrange to perform scREALTIME----
collect_expression_all = NULL
for (i in seq(1:nrow(param_sets))){
  print(i)
  collect_expression_paramset = NULL
  
  for (j in  c("AP1","NFkB","NFkB|p38","NFkB|IRF","IRF")){
    print(j)
    solve = readRDS(paste0("./generated_scTFinputs/output_test_5GRS_paramsets/output_test_5GRS_paramset",i,".rds"))
    #label each gene by GRS_n_Kd_kdeg
    # solve$genename = paste0(solve$GRS,"_", solve$params[["n"]], "_",
    #                         solve$params[["ktA"]], "_", signif(solve$params[["k_deg"]],2))
    tmp = solve[solve$GRS==j,] #separate out each GRS
    colnames(tmp)[2] = paste0(j,"_n", param_sets$n[i], "_ktA",
                              param_sets$ktA[i], "_kdeg", log(2)/param_sets$k_deg[i])
    
    if(is.null(collect_expression_paramset)){
      collect_expression_paramset = rbind(collect_expression_paramset, tmp[,1:(ncol(tmp)-2)]) #remove paramset and GRS
    }else{
      collect_expression_paramset = cbind(collect_expression_paramset, data.frame(tmp[,2,drop=F])) # add mRNA col
    }
  }
  collect_expression_paramset = collect_expression_paramset[c(1,3,4,2,5:ncol(collect_expression_paramset))]
  
  if(is.null(collect_expression_all)){
    collect_expression_all = rbind(collect_expression_all, collect_expression_paramset)
  }else{
    collect_expression_all = cbind(collect_expression_all, collect_expression_paramset[,-c(1:3)])
  }
  
}
# saveRDS(collect_expression_all, paste0("./generated_scTFinputs/output_test_5GRS_all40paramsets_LPS.rds"))
# saveRDS(collect_expression_all, paste0("./generated_scTFinputs/output_test_5GRS_all70paramsets_allstims.rds"))

#################################----
#rearrange to perform scREALTiME and add 10% transcriptional noise at the same time-----
collect_expression_all = NULL
for (i in seq(1:nrow(param_sets))){
  print(i)
  collect_expression_paramset = NULL
  
  for (j in  c("AP1","NFkB","NFkB|p38","NFkB|IRF","IRF")){
    print(j)
    solve = readRDS(paste0("./generated_scTFinputs/output_test_5GRS_paramset",i,"_LPS.rds"))
    #label each gene by GRS_n_Kd_kdeg
    # solve$genename = paste0(solve$GRS,"_", solve$params[["n"]], "_",
    #                         solve$params[["ktA"]], "_", signif(solve$params[["k_deg"]],2))
    tmp = solve[solve$GRS==j,] #separate out each GRS
    colnames(tmp)[2] = paste0(j,"_n", param_sets$n[i], "_ktA",
                              param_sets$ktA[i], "_kdeg", log(2)/param_sets$k_deg[i])
    
    if(is.null(collect_expression_paramset)){
      collect_expression_paramset = rbind(collect_expression_paramset, tmp[,1:(ncol(tmp)-2)]) #remove paramset and GRS
    }else{
      collect_expression_paramset = cbind(collect_expression_paramset, data.frame(tmp[,2,drop=F])) # add mRNA col
    }
  }
  collect_expression_paramset = collect_expression_paramset[c(1,3,4,2,5:ncol(collect_expression_paramset))]
  
  if(is.null(collect_expression_all)){
    collect_expression_all = rbind(collect_expression_all, collect_expression_paramset)
  }else{
    collect_expression_all = cbind(collect_expression_all, collect_expression_paramset[,-c(1:3)])
  }
  
}

jitter.custom <- function(x) jitter(x, factor = 5, amount = 0)

collect_expression_jitter = cbind(collect_expression_all[,c(1:3)],
                                  data.frame(lapply(collect_expression_all[,-c(1:3)], FUN = jitter.custom)))
ggplot(collect_expression_jitter, aes(time,`NFkB_n1_ktA0.48_kdeg30`, group = cell)) +
  # geom_vline(xintercept = c(0,0.25,1,3,8), linetype="dotted")+ #theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  # geom_vline(xintercept = c(0,0.5,1,3,5,8), linetype="dotted")+
  geom_line(aes(group = as.factor(cell), color = stimulus), alpha = 0.5)+
  theme_bw(base_size = 14)+theme(legend.position = "none")
# saveRDS(collect_expression_jitter, paste0("./generated_scTFinputs/output_txnnoisejitter_5GRS_all40paramsets_LPS.rds"))
# saveRDS(collect_expression_jitter, paste0("./generated_scTFinputs/output_txnnoise10p.jitter_5GRS_all40paramsets_LPS.rds"))


#make Seurat object----
library(scREALTIME)
if(1){
  # collect_expression_all = readRDS("./generated_scTFinputs/output_txnnoisejitter_5GRS_all40paramsets_LPS.rds")
  # collect_expression_all = readRDS("./generated_scTFinputs/output_txnnoise10p.jitter_5GRS_all40paramsets_LPS.rds")
  # collect_expression_all = readRDS("./generated_scTFinputs/output_test_5GRS_all40paramsets_LPS.rds")
  
  timepoints = unique(collect_expression_all$time)
  num_used_timepoints =5 # Use 5 evenly spaced timepoints
  # 480 min = 8 hours, 0 hr time point, integrated over 300 timepoints
  
  select_timepoints = seq(1,length(timepoints),length.out = num_used_timepoints)
  select_timepoints = seq(2,length(timepoints),length.out = num_used_timepoints)#avoiding 0 timept for simulated data
  select_timepoints = timepoints[select_timepoints]
  
  #handpick select timepoints
  # select_timepoints = timepoints[c(2,20,39,113,300)] #0,0.5hr,1hr,3hr, 8hr
  # select_timepoints = timepoints[c(2,188,263,281,300)] #0, 5hr,7hr, 7.5hr, 8hr
  
  print(timepoints)
  print(select_timepoints)
  
  # Show PCA trajectories and true trajectories (for each gene)
  num_cells_for_lineplots = 50 #50 simulated cells
  numgenes = 200 #(40paramsets * 5GRSs)
  
  pca = prcomp(collect_expression_all[1:(length(timepoints) * num_cells_for_lineplots),4:(numgenes+3)], center = F,scale = FALSE) # 3 cells
  inds = seq(1,(length(timepoints) * num_cells_for_lineplots), length.out = num_used_timepoints * num_cells_for_lineplots)
  
  # Cast to Seurat object
  
  data_for_Seurat = as.data.frame(t(collect_expression_all[collect_expression_all$time %in% select_timepoints,4:(numgenes+3)]))
  meta = collect_expression_all[collect_expression_all$time %in% select_timepoints,1:3]
  rownames(meta)= paste0(meta$label,"_",meta$cell,"_",round(meta$time,2))
  colnames(data_for_Seurat) = rownames(meta)
  # transpose
  # store in RNA assay slot directly -- no normalization
  require('Seurat')
  macro_LPS= CreateSeuratObject(data_for_Seurat, meta.data = meta, project = 'simulations')
  macro_LPS@meta.data$stimulus = sapply(strsplit(rownames(macro_LPS@meta.data),"_"), `[`, 1)
  macro_LPS@meta.data$cell = paste0("cell_",sapply(strsplit(rownames(macro_LPS@meta.data),"_"), `[`, 3))
  macro_LPS@meta.data$timept = paste(sapply(strsplit(rownames(macro_LPS@meta.data),"_"), `[`, 4), 'min', sep = '')
  macro_LPS@meta.data$timept_num = as.numeric(sapply(strsplit(rownames(macro_LPS@meta.data),"_"), `[`, 4))
  head(macro_LPS@meta.data)
  # Data viz
  # VlnPlot(macro_LPS, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
  # FeatureScatter(macro_LPS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  #Finding variable features (note: normalization skipped)
  macro_LPS <- FindVariableFeatures(macro_LPS, selection.method = "vst", nfeatures = 200)
  
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(macro_LPS), 10)
  top10
  
  # plot variable features with and without labels
  plot1 <- VariableFeaturePlot(macro_LPS)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  plot2
  
  if(0){
    #plot variable features over timepts
    collect_variablefeatures = NULL
    for (t in unique(macro_LPS$timept_num)){
      print(t)
      macro_LPS.subset = subset(macro_LPS, subset = timept_num ==t)
      macro_LPS.subset = FindVariableFeatures(macro_LPS.subset, selection.method = "vst", nfeatures = 200)
      tmp = macro_LPS.subset$RNA@meta.features
      collect_variablefeatures = cbind(collect_variablefeatures, tmp$vst.variance)
    }
    collect_variablefeatures = as.data.frame(collect_variablefeatures)
    colnames(collect_variablefeatures) = unique(macro_LPS$timept_num)
    rownames(collect_variablefeatures) = rownames(macro_LPS$RNA@meta.features)
    collect_variablefeatures = cbind(gene = rownames(collect_variablefeatures), collect_variablefeatures)
    collect_variablefeatures.m = melt(collect_variablefeatures)
    collect_variablefeatures.m$GRS = gsub("-..*","",collect_variablefeatures.m$gene)
    ggplot(collect_variablefeatures.m, aes(as.numeric(variable), log10(value)))+
      geom_path(aes(group = gene, color = GRS), alpha = 0.25)+ ylab("log10(variance)")+xlab("time")+
      geom_point(aes(color = GRS),size=0.5)+theme_classic(base_size = 14)+facet_wrap(~GRS,ncol = 3)
  }
  
  # scaling (before PCA)
  macro_LPS <- FindVariableFeatures(macro_LPS)
  all.genes <- rownames(macro_LPS)
  macro_LPS <- ScaleData(macro_LPS)
  
  #PCA
  macro_LPS[["RNA"]]@scale.data = as.matrix(macro_LPS[["RNA"]]@data)
  macro_LPS <- RunPCA(macro_LPS,assay = "RNA", scale =F, center = F)
  print(macro_LPS[["pca"]], dims = 1:5, nfeatures = 5)
  VizDimLoadings(macro_LPS, dims = 1:2, reduction = "pca")
  DimPlot(macro_LPS, reduction = "pca", group.by = "timept")#+theme(legend.position = "None")
  DimHeatmap(macro_LPS, dims = 1:5, cells = 50, balanced = TRUE)
  
  #plot PCA with lines connecting timepts
  data = macro_LPS@meta.data
  pca = data.frame(macro_LPS@reductions$pca@cell.embeddings)
  data$PC1 = pca$PC_1[match(rownames(data), rownames(macro_LPS@reductions$pca))]
  data$PC2 = pca$PC_2[match(rownames(data), rownames(macro_LPS@reductions$pca))]
  head(data)
  ggplot(data, aes(PC1,PC2, group=cell))+geom_path(aes(group = cell), alpha=0.25)+
    geom_point(aes(color = timept),size=1.5)+theme_classic(base_size = 16)+theme(legend.position = "none")
  
  macro_LPS@active.assay
  # saveRDS(macro_LPS, paste0("./generated_scTFinputs/LPS_simulatedwGeneratedscTF_",num_used_timepoints,"timepts.rds"))
  # saveRDS(macro_LPS, paste0("./generated_scTFinputs/LPS_simulatedwGeneratedscTF.txnnoisejitter_",num_used_timepoints,"timepts.rds"))
  # saveRDS(macro_LPS, paste0("./generated_scTFinputs/LPS_simulatedwGeneratedscTF.txnnoise10p.jitter_",num_used_timepoints,"timepts.rds"))
  
  #5timepts as measured
  # saveRDS(macro_LPS, paste0("./generated_scTFinputs/LPS_simulatedwGeneratedscTF_",num_used_timepoints,"timepts.selected.rds"))
  # saveRDS(macro_LPS, paste0("./generated_scTFinputs/LPS_simulatedwGeneratedscTF.txnnoisejitter_",num_used_timepoints,"timepts.selected.rds"))
  # saveRDS(macro_LPS, paste0("./generated_scTFinputs/LPS_simulatedwGeneratedscTF.txnnoise10p.jitter_",num_used_timepoints,"timepts.selected.rds"))
  
}

# run scREALTIME
library(scREALTIME)
if (1){
  
  num_used_timepoints = 5
  macro_LPS = readRDS(paste0("./generated_scTFinputs/LPS_simulatedwGeneratedscTF_",num_used_timepoints,"timepts.rds"))
  
  #5timepts selected
  # macro_LPS = readRDS(paste0("./generated_scTFinputs/LPS_simulatedwGeneratedscTF_",num_used_timepoints,"timepts.selected.rds"))
  
  select_timepoints = unique(macro_LPS$time)
  
  ## Sims
  n="LPS"
  select_timepoints = round(select_timepoints,2)
  PCAPlot(macro_LPS, dims = c(1,2), cells = as.vector(rownames(metadata)), group.by = 'timept', shuffle = TRUE, label = TRUE, label.box = TRUE) + ggtitle('PCA with PC 2 and 3')
  
  num_archetypes = 20 #if no file label, then default was 10 archetypes
  interpl = "spline.mono"
  
  set.seed(1)
  reconst = getTrajectory4(macro = macro_LPS, metadata = macro_LPS@meta.data, num_archetypes = num_archetypes,data = "RNA",
                           timepoints = select_timepoints, num_trajectories = 100, num_sim_pts = 100,
                           reduction = 'pca', stimulus = n, consensus_measure = 'mean',
                           interpolant = interpl, prob_method = 'distance', distance_metric = 'euclidean' ,
                           varFilter = F, exp_prob = 1) 
  
  saveRDS(reconst, paste0("./generated_scTFinputs/LPS_",interpl,"_uncentered_",num_used_timepoints,"timepts.used_",num_archetypes,"archetypes.rds"))
  
  #5timepts selected
  # saveRDS(reconst, paste0("./generated_scTFinputs/LPS_",interpl,"_uncentered_",num_used_timepoints,"timepts.selected.used_",num_archetypes,"archetypes.rds"))
  
}


#####################################################################
#Figure 2: Evaluation of scGETs imputation method against ground truth
#####################################################################

#test the probability matrix step without k-means---------------------------
# install.packages("ggalluvial")
library(ggalluvial)

metrics = matrix(nrow = 50) #50 ground truth trajectories

for(num_timepts in c(5,10,20,50)){
  for(num_comps in c(2,5,10,20,50)){ #only run inner loop for 5.selected
    if (1){
      num_used_timepoints = 5 #num_timepts
      print(num_used_timepoints)
      num_comps = 50
      print(num_comps)
      
      #5timepts selected
      # macro_LPS = readRDS(paste0("./generated_scTFinputs/LPS_simulatedwGeneratedscTF_",num_used_timepoints,"timepts.selected.rds"))
      # macro_LPS = readRDS(paste0("./generated_scTFinputs/LPS_simulatedwGeneratedscTF.txnnoisejitter_",num_used_timepoints,"timepts.selected.rds"))
      # macro_LPS = readRDS(paste0("./generated_scTFinputs/LPS_simulatedwGeneratedscTF.txnnoise10p.jitter_",num_used_timepoints,"timepts.selected.rds"))
      
      # macro_LPS = readRDS(paste0("./generated_scTFinputs/LPS_simulatedwGeneratedscTF_",num_used_timepoints,"timepts.rds"))
      # macro_LPS = readRDS(paste0("./generated_scTFinputs/LPS_simulatedwGeneratedscTF.txnnoisejitter_",num_used_timepoints,"timepts.rds"))
      # macro_LPS = readRDS(paste0("./generated_scTFinputs/LPS_simulatedwGeneratedscTF.txnnoise10p.jitter_",num_used_timepoints,"timepts.rds"))
      DimPlot(macro_LPS, reduction = "pca", group.by = "timept")+theme(legend.position = "None")
      
      select_timepoints = unique(macro_LPS$timept_num)
      
      num_archetypes = 50 #if no file label, then default was 10 archetypes
      interpl = "spline"
      n = "LPS"
      
      set.seed(1)
      macro = macro_LPS; metadata = macro_LPS@meta.data; num_archetypes = num_archetypes;data = "RNA"
      timepoints = select_timepoints; num_trajectories = 500; num_sim_pts = 100;
      reduction = 'pca'; stimulus = n; consensus_measure = 'mean'; 
      interpolant = interpl; prob_method = 'distance'; distance_metric = 'euclidean' 
      varFilter = F;exp_prob = 1
      pc_comps=num_comps;
      
      if (1){
        cells_by_timept <- list()
        for(i in timepoints){
          index = paste("time_",i, "hr", sep = "")
          query = paste(i, "hr", sep='')
          #cells_by_timept[[index]] <- rownames(metadata)[metadata$timept == query]
          cells_by_timept[[index]] <- rownames(metadata)[metadata$timept_num == i]
        }
        
        if(data == 'RNA'){
          RNA <- as.data.frame(macro@assays$RNA@data)
        }else if(data == 'ISNorm'){
          RNA <- as.data.frame(macro@assays$ISnorm@data)
        }
        
        RNA <- t(RNA)
        RNA <- RNA[rownames(RNA) %in% rownames(metadata),]
        
        clusterings <- list()
        zero_var_inds <- list()
        
        for(i in timepoints){
          index = paste("time_", i, "hr", sep = "")
          data <- RNA[rownames(RNA) %in% cells_by_timept[[index]],]
          zero_var_inds[[index]] <- resample::colVars(data) == 0
          print(zero_var_inds)
          print(tail(data[,! zero_var_inds[[index]]]))
          print(dim(data[,! zero_var_inds[[index]]]))
          
          # if(varFilter){
          #   clusterings[[index]] <- kmeans(data[,!zero_var_inds[[index]]], centers = num_archetypes, iter.max = 50)
          # }else{
          #   clusterings[[index]] <- kmeans(data, centers = num_archetypes, iter.max = 50)
          # }
          
          if(1){ #no kmeans for simulations
            clusterings[[index]] <- list(cluster=(seq(1:num_archetypes))) #keep original
            clusterings[[index]] <- list(cluster=(sample(1:num_archetypes, num_archetypes, F))) #shuffle connections
            names(clusterings[[index]]$cluster) = c(cells_by_timept[[index]])
          }
          
          if(i == timepoints[1]){
            cluster_counts <- as.data.frame(table(clusterings[[index]]$cluster))
          }else{
            cluster_counts <- cbind(cluster_counts, table(clusterings[[index]]$cluster))
          }
        }
        
        ## Clerical in order to make cluster counts data frame more clean
        col_omit <- c(1:length(timepoints))*2 - 1
        cluster_counts <- cluster_counts[, -col_omit]
        rownames(cluster_counts) <- c(1:num_archetypes)
        rownames(cluster_counts) <- paste("bin", rownames(cluster_counts), "")
        colnames(cluster_counts) <- timepoints
        colnames(cluster_counts) <- paste("time_", colnames(cluster_counts), sep="")
        
        DimPlot(macro, reduction = "pca", group.by = "timept")
        if(toupper(reduction) == 'PCA'){
          # pca = prcomp(RNA, center = F, scale = F, rank. = 50)
          # pcscores = pca$x
          # loadings = pca$rotation
          pcscores = macro[['pca']]@cell.embeddings[,1:pc_comps]
          #print(pcscores)
        }else if(toupper(reduction) == 'NMF'){
          pcscores = macro[['NMF']]@cell.embeddings[,1:pc_comps]
        }else if(toupper(reduction) == 'ICA'){
          pcscores = macro[['ica']]@cell.embeddings[,1:pc_comps]
        }
        
        
        #plot for sanity
        require('factoextra')
        require('ggfortify')
        
        #fviz_eig(pca)
        #autoplot(pca, data = df_for_pca, loadings = FALSE, colour = 'metadata.timept')
        
        cluster_densities <- as.data.frame(prop.table(as.matrix(cluster_counts), 2))
        cluster_densities <- cluster_densities^(exp_prob)
        
        cell_cluster_df <- NULL
        for(i in timepoints){
          index = paste("time_", i, "hr", sep = "")
          cell_cluster_df <- c(cell_cluster_df, clusterings[[index]]$cluster)
        }
        names = names(cell_cluster_df)
        cell_cluster_df <- as.data.frame(cell_cluster_df)
        rownames(cell_cluster_df)= names
        
        pcscores <- as.data.frame(pcscores)
        pcscores_stim <- pcscores[rownames(pcscores) %in% rownames(metadata),]
        pcscores_stim$timept <- metadata$timept_num[match(rownames(pcscores_stim), rownames(metadata))]
        pcscores_stim$bin <- cell_cluster_df$cell_cluster_df[match(rownames(pcscores_stim), rownames(cell_cluster_df))]
        pcscores_stim$timebin_tag <- paste(pcscores_stim$timept, pcscores_stim$bin, sep = "_")
        
        if(consensus_measure == 'mean'){
          mean_pcscores <- aggregate(pcscores_stim[, 1:(ncol(pcscores_stim)-3)], list(pcscores_stim$timebin_tag), mean)
        }else if(consensus_measure == 'median'){
          mean_pcscores <- aggregate(pcscores_stim[, 1:(ncol(pcscores_stim)-3)], list(pcscores_stim$timebin_tag), median)
        }
        mean_pcscores$timept <- sapply(strsplit(mean_pcscores$Group.1, split = "_"), `[`, 1)
        mean_pcscores$bin <- sapply(strsplit(mean_pcscores$Group.1, split = "_"), `[`, 2)
        col_orders = c(1,(ncol(mean_pcscores)), (ncol(mean_pcscores)-1), 2:(ncol(mean_pcscores)-2))
        mean_pcscores <- mean_pcscores[,col_orders]
        
        mean_pcscores = mean_pcscores[order(as.factor(as.numeric(mean_pcscores$timept))),]
        ggplot(mean_pcscores, aes(PC_1,PC_2, group=bin))+geom_path(aes(group = bin), alpha=0.25)+
          geom_point(aes(color = as.factor(as.numeric(timept))),size=1.5)+theme_bw(base_size = 14)+theme(legend.position = "none")
        
        #View(mean_pcscores)
        #print(str(mean_pcscores))
        
        # Generating walk probabilities matrix
        
        # data = macro_LPS@meta.data
        # pca = data.frame(macro_LPS@reductions$pca@cell.embeddings)
        # data = cbind(data, pca[match(rownames(data), rownames(macro_LPS@reductions$pca)),])
        # head(data)
        # ggplot(data, aes(PC_1,PC_2, group=cell))+geom_path(aes(group = cell), alpha=0.25)+
        #   geom_point(aes(color = timept),size=1.5)+theme_bw(base_size = 12)+theme(legend.position = "none")
        
        
        distances = as.matrix(dist(mean_pcscores[,4:ncol(mean_pcscores)]))
        rownames(distances) = mean_pcscores$Group.1
        colnames(distances) = mean_pcscores$Group.1
        # distances = distances[order((as.numeric(gsub("..*_","",rownames(distances))))),
        #                       order((as.numeric(gsub("..*_","",colnames(distances)))))]
        # p1=pheatmap(log10(distances[grepl("^1.6",rownames(distances)), grepl("^30",colnames(distances))]+1), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
        # p2=pheatmap(log10(distances[grepl("^30",rownames(distances)), grepl("^61",colnames(distances))]+1), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
        # p3=pheatmap(log10(distances[grepl("^61",rownames(distances)), grepl("^179",colnames(distances))]+1), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
        # p4=pheatmap(log10(distances[grepl("^179",rownames(distances)), grepl("^480",colnames(distances))]+1), cluster_rows = F, cluster_cols = F, show_rownames = F, show_colnames = F)
        
        # Need to convert this to a transitions probs matrix at each step
        
        # Need to normalize so larger distances have lower probability
        
        
        
        walks <- list()
        # walk_probs <- matrix(nrow = num_trajectories, ncol = 2)
        # walk_probs[,1] = c(1:num_trajectories)
        set.seed(123)
        if(prob_method == 'distance'){
          print("Generating random walks based on distance between clusters")
          for(i in 1:num_trajectories){
            # cur_prob = 1
            path <- c()
            for(j in 1:length(timepoints)){
              if(j==1){
                next_value = sample(c(1:num_archetypes), size = 1, prob = cluster_densities[,j])
              }
              else{
                prev = path[length(path)]
                row = (num_archetypes)*(j-2) + prev
                cols = (num_archetypes*(j-1) + 1):(num_archetypes*j)
                probs = distances[row, cols]
                #probs = exp(sum(probs)-probs)
                #probs = sum(probs)-probs
                probs = (1/probs)^exp_prob
                #print(probs/sum(probs))
                next_value = sample(as.numeric(mean_pcscores$bin[1:num_archetypes]), size = 1, prob = probs) #as.numeric(mean_pcscores$bin[1:50]) #as.numeric(sub(".*?_","",names(probs)[probs==max(probs)]))#
                # ggplot(data.frame(archetype =sample(as.numeric(mean_pcscores$bin[1:num_archetypes]), 
                #                                     size = 5000, replace = T, prob = probs)), aes(archetype))+
                #   geom_histogram(bins = 50)
              }
              # cur_prob = cur_prob * cluster_densities[next_value,j]
              path <- c(path, next_value)
            }
            walks[[i]] <- path
            # walk_probs[i,2] <- cur_prob
          }
        }else if(prob_method == 'density'){
          print("Generating random walks based on cluster densities")
          walks <- list()
          walk_probs <- matrix(nrow = num_trajectories, ncol = 2)
          walk_probs[,1] = c(1:num_trajectories)
          for(i in 1:num_trajectories){
            cur_prob = 1
            path <- c()
            for(j in 1:length(timepoints)){
              next_value = sample(c(1:num_archetypes), size = 1, prob = cluster_densities[,j])
              cur_prob = cur_prob * cluster_densities[next_value,j]
              path <- c(path, next_value)
            }
            walks[[i]] <- path
            walk_probs[i,2] <- cur_prob
          }
        }else if(prob_method == 'hybrid'){
          for(i in 1:num_trajectories){
            cur_prob = 1
            path <- c()
            for(j in 1:length(timepoints)){
              if(j==1){
                next_value = sample(c(1:num_archetypes), size = 1, prob = cluster_densities[,j])
              }
              else{
                prev = path[length(path)]
                row = (num_archetypes)*(j-2) + prev
                cols = (num_archetypes*(j-1) + 1):(num_archetypes*j)
                probs = distances[row, cols]
                #probs = exp(sum(probs)-probs)
                #probs = sum(probs)-probs
                probs = (1/probs)^exp_prob
                #print(probs/sum(probs))
                next_value = sample(c(1:num_archetypes), size = 1, prob = probs)
              }
              cur_prob = cur_prob * cluster_densities[next_value,j]
              path <- c(path, next_value)
            }
            walks[[i]] <- path
            walk_probs[i,2] <- cur_prob
          }
        }else{
          print('Invalid method for random walk probabilities')
        }
        
        # Walks now hold bins for random walks
        # Checking number of unique walks
        num_unique_traj = length(unique(sapply( walks, paste0, collapse="")))
        print(paste(length(unique(sapply( walks, paste0, collapse=""))), "unique walks/trajectories"))
        
        # Make table of unique walks for us to keep track of
        walk_frequencies <- as.data.frame(table(sapply( walks, paste0, collapse="")))
        colnames(walk_frequencies)[1] = 'Walk'
        
        #plot histograms
        walk.df = data.frame(t(as.data.frame(walks)))
        colnames(walk.df) = select_timepoints
        # ggplot(walk.df, aes(X1))+ geom_histogram(color="black", fill="gray", bins = 50)+
        #   theme_bw(base_size = 14)
        walk.df.m = melt(walk.df)
        ggplot(walk.df.m, aes(value))+ geom_histogram( fill="gray", color = "black", bins = 50)+
          theme_bw(base_size = 10)+facet_grid(~variable)+xlab("archetype")
        
        walk.cast = dcast(walk.df.m, value~variable)
        colnames(walk.cast)[1] = "bin"
        
        
        spline_pts <- list()
        for(pc in 1:(ncol(mean_pcscores)-3)){
          id = paste("pc", pc, sep="_")
          spline_pts[[id]] <- matrix(nrow = num_trajectories, ncol = length(timepoints))
          #View(spline_pts)
          for(tr in 1:num_trajectories){
            path = walks[[tr]]
            for(j in 1:length(timepoints)){
              time = timepoints[j]
              bin = path[j]
              spline_pts[[id]][tr, j] <- mean_pcscores[round(as.numeric(mean_pcscores$timept),3) == round(time,3) & mean_pcscores$bin == bin, pc+3]
            }
          }
        }
        
        
        spline_pts_PC1 = (spline_pts$pc_1)
        colnames(spline_pts_PC1) = timepoints
        spline_pts_PC1.m = melt(spline_pts_PC1)
        spline_pts_PC1.m$cell = paste0("archetype_",spline_pts_PC1.m$Var1,"_",spline_pts_PC1.m$Var2)
        colnames(spline_pts_PC1.m)[1:3] = c("bin","timept", "pc_1")
        
        for(pc in 2:(ncol(mean_pcscores)-3)){
          id = paste("pc", pc, sep="_")
          print(id)
          
          spline_pts_PCx = spline_pts[[id]]
          colnames(spline_pts_PCx) = timepoints
          spline_pts_PCx.m = melt(spline_pts_PCx)
          spline_pts_PCx.m$cell = paste0("archetype_",spline_pts_PCx.m$Var1,"_",spline_pts_PCx.m$Var2)
          colnames(spline_pts_PCx.m)[1:2] = c("bin","timept")
          
          spline_pts_PC1.m = cbind(spline_pts_PC1.m, spline_pts_PCx.m$value[match(spline_pts_PC1.m$cell, spline_pts_PCx.m$cell)] )
          colnames(spline_pts_PC1.m)[ncol(spline_pts_PC1.m)] = id
        }
        ggplot(spline_pts_PC1.m, aes(pc_1,pc_2, group=bin))+geom_path(aes(group = bin), alpha=0.25)+
          geom_point(aes(color = as.factor(timept)),size=1.5)+theme_bw(base_size = 14)+theme(legend.position = "none")
        
        #plot PCs vs time
        # ggplot(spline_pts_PC1.m, aes(timept, pc_1*-1))+geom_path(aes(group = bin), alpha=0.25)+
        #   geom_point(aes(color = as.factor(timept)),size=1.5)+theme_bw(base_size = 14)+theme(legend.position = "none")
        # ggplot(spline_pts_PC1.m, aes(timept, pc_2))+geom_path(aes(group = bin), alpha=0.25)+
        #   geom_point(aes(color = as.factor(timept)),size=1.5)+theme_bw(base_size = 14)+theme(legend.position = "none")
        # ggplot(spline_pts_PC1.m, aes(timept, pc_3))+geom_path(aes(group = bin), alpha=0.25)+
        #   geom_point(aes(color = as.factor(timept)),size=1.5)+theme_bw(base_size = 14)+theme(legend.position = "none")
        # 
        
        if(1){ #spline fitting
          
          
          # Now, we get number of unique paths * number of trajectories * 50 PCs
          num_splines =dim(walk_frequencies)[1] * pc_comps
          print(paste("Number of splines: ", dim(walk_frequencies)[1] * pc_comps))
          # dup_rows = duplicated(spline_pts[[1]])
          # walk_probs = walk_probs[!dup_rows, ]
          # 
          # for(i in 1:length(spline_pts)){
          #   spline_pts[[i]] <- spline_pts[[i]][!dup_rows,]
          #   
          #   #spline_pts[[i]] <- unique(spline_pts[[i]])
          # }
          
          ## Iterate over all PCs
          #num_sim_pts = 100
          sim_times = seq(min(timepoints),max(timepoints),length.out = num_sim_pts)
          for(i in timepoints){ # TIMEPOINTS
            if(!(i %in% sim_times)){
              sim_times <- c(sim_times, i)
              num_sim_pts = num_sim_pts + 1
            }
          }
          sim_times <- sort(sim_times)
          
          simulated = matrix(ncol = (ncol(mean_pcscores)-3), nrow = num_sim_pts*sum(walk_frequencies$Freq))
          simulation_mat_time = rep(sim_times, sum(walk_frequencies$Freq))
          
          for(i in 1:(ncol(mean_pcscores)-3)){ # Loop over PCs
            preds = c()
            for(walk in 1:nrow(spline_pts[[i]])){
              if(interpolant == 'spline'){
                spline <- smooth.spline(timepoints, spline_pts[[i]][walk,], cv = T, all.knots = T) # TIMEPOINTS
                preds <- c(preds, predict(spline, sim_times)$y)
              } else if(interpolant == 'spline.df3'){
                spline <- smooth.spline(timepoints, spline_pts[[i]][walk,], df=3, all.knots = T) # TIMEPOINTS
                preds <- c(preds, predict(spline, sim_times)$y)
              } else if(interpolant == 'spline.df2'){
                spline <- smooth.spline(timepoints, spline_pts[[i]][walk,], df=2, all.knots = T) # TIMEPOINTS
                preds <- c(preds, predict(spline, sim_times)$y)
              } else if(interpolant == 'spline.df1'){
                spline <- smooth.spline(timepoints, spline_pts[[i]][walk,], df=1, all.knots = T) # TIMEPOINTS
                preds <- c(preds, predict(spline, sim_times)$y)
              } else if(interpolant == 'spline.df5'){
                spline <- smooth.spline(timepoints, spline_pts[[i]][walk,], df=5, all.knots = T) # TIMEPOINTS
                preds <- c(preds, predict(spline, sim_times)$y)
              } else if(interpolant == 'spline.mono'){
                fun = splinefun(timepoints, spline_pts[[i]][walk,], method = "monoH.FC")
                preds = c(preds, fun(sim_times))
              } else if(interpolant == 'spline.periodic'){
                fun = splinefun(timepoints, spline_pts[[i]][walk,], method = "periodic")
                preds = c(preds, fun(sim_times))
              } else if(interpolant == 'spline.natural'){
                fun = splinefun(timepoints, spline_pts[[i]][walk,], method = "natural")
                preds = c(preds, fun(sim_times))
              } else if(interpolant == 'linear'){
                fun = approxfun(timepoints, spline_pts[[i]][walk,], rule = 2)
                preds = c(preds, fun(sim_times))
              } else if(interpolant == 'loess'){
                loe = loess(spline_pts[[i]][walk,] ~ timepoints, span = 0.75, se = F)
                preds = c(preds, predict(loe, sim_times))
              }else {
                print('Unspecified interpolant')
              }
            }
            simulated[,i] = preds
          }
          
        }
        
        #plot spline fit per PC
        if(0){
          data = macro_LPS@meta.data
          pca = data.frame(macro_LPS@reductions$pca@cell.embeddings)
          data= cbind(data, pca[ match(rownames(data), rownames(macro_LPS@reductions$pca)),])
          
          simulated.annot = data.frame(simulated)
          simulated.annot$time <- simulation_mat_time
          simulated.annot$path <- rep(1:sum(walk_frequencies$Freq), each = num_sim_pts)
          
          ggplot(simulated.annot[grepl("^1$", simulated.annot$path),], aes(time, X1))+geom_path()+
            geom_point(data[grepl("cell_79",data$cell),], mapping = aes(timept_num,PC_1, color = as.factor(time)),size=3)+
            theme_bw(base_size = 12)+ylab("PC1")
          
          ggplot(simulated.annot[grepl("^1$", simulated.annot$path),], aes(time, X2))+geom_path()+
            geom_point(data[grepl("",data$cell),], mapping = aes(timept_num,PC_2, color = as.factor(time)),size=1)+
            theme_bw(base_size = 12)+ylab("PC2")
          
          ggplot(simulated.annot[grepl("^1$", simulated.annot$path),], aes(X1, X2))+geom_path()+
            geom_point(data[grepl("",data$cell),], mapping = aes(PC_1,PC_2, color = as.factor(time)),size=1)+
            theme_bw(base_size = 12)+xlab("PC1")+ylab("PC2")
          
          ggplot(simulated.annot[grepl("", simulated.annot$path),], aes(X1, X2))+geom_path(aes(group=path),alpha=0.05)+
            geom_point(data[grepl("",data$cell),], mapping = aes(PC_1,PC_2, color = as.factor(time)),size=1)+
            theme_bw(base_size = 12)+xlab("PC1")+ylab("PC2")
        }
        
        # Now let's cast back to gene expression space
        loadings <- macro[[reduction]]@feature.loadings[,1:pc_comps]
        
        if(toupper(reduction) == 'NMF'){
          loadings = macro[['NMF']]@feature.loadings[,1:pc_comps]
        }
        if(toupper(reduction) == 'ICA'){
          loadings = macro[['ica']]@feature.loadings[,1:pc_comps]
        }
        #dim(loadings)
        reconstructed_pc <- t(loadings %*% t(simulated))
        reconstructed_pc <- as.data.frame(reconstructed_pc)
        reconstructed_pc$time <- simulation_mat_time
        reconstructed_pc$path <- rep(1:sum(walk_frequencies$Freq), each = num_sim_pts)
        # 
        toRet = list()
        toRet[['reconstructed_trajectories']] = reconstructed_pc
        toRet[['cluster_densities']] = cluster_densities
        toRet[['metadata']] = metadata
        toRet[['number_unique_trajectories']] = num_unique_traj
        toRet[['number_of_splines']] = num_splines
        # toRet[['probability_of_walks']] = walk_probs
        
        call = list()
        call[['Seurat_obj']] = macro
        call[['metadata_name']] = metadata
        call[['num_archetypes']] = num_archetypes
        call[['timepoints']] = timepoints
        call[['num_trajectories']] = num_trajectories
        call[['num_sim_pts']] = num_sim_pts
        call[['reduction']] = reduction
        call[['stimulus']] = stimulus
        call[['consensus_measure']] = consensus_measure
        toRet[['call']] = call
        
        # return(toRet)
      }
      
      #5selectedtimepts
      # saveRDS(toRet, paste0("./generated_scTFinputs/LPS_",interpl,"_uncentered_",num_used_timepoints,"timepts.used.selected_",num_archetypes,"archetypes_nokmeans.rds"))
      # saveRDS(toRet, paste0("./generated_scTFinputs/LPS.txnnoise_",interpl,"_uncentered_",num_used_timepoints,"timepts.used.selected_",num_archetypes,"archetypes_nokmeans.rds"))
      # saveRDS(toRet, paste0("./generated_scTFinputs/LPS.txnnoise10p_",interpl,"_uncentered_",num_used_timepoints,"timepts.used.selected_",num_archetypes,"archetypes_nokmeans.rds"))
      
      # saveRDS(toRet, paste0("./generated_scTFinputs/LPS_",interpl,"_uncentered_",num_used_timepoints,"timepts.used_",num_archetypes,"archetypes_nokmeans.rds"))
      # saveRDS(toRet, paste0("./generated_scTFinputs/LPS.txnnoise_",interpl,"_uncentered_",num_used_timepoints,"timepts.used_",num_archetypes,"archetypes_nokmeans.rds"))
      # saveRDS(toRet, paste0("./generated_scTFinputs/LPS.txnnoise10p_",interpl,"_uncentered_",num_used_timepoints,"timepts.used_",num_archetypes,"archetypes_nokmeans.rds"))
      
      #only 50trajectories, 500traj above
      # saveRDS(toRet, paste0("./generated_scTFinputs/LPS_",interpl,"_uncentered_",num_used_timepoints,"timepts.used.selected_",num_archetypes,"archetypes_nokmeans_50traj.rds"))
      # saveRDS(toRet, paste0("./generated_scTFinputs/LPS_",interpl,"_uncentered_",num_used_timepoints,"timepts.used_",num_archetypes,"archetypes_nokmeans_50traj.rds"))
      
      #using exp_prob=xx
      # saveRDS(toRet, paste0("./generated_scTFinputs/LPS_",interpl,"_uncentered_",num_used_timepoints,"timepts.used_",num_archetypes,"archetypes_nokmeans_exp",exp_prob,".rds"))
      
      #using pc_comps==xx
      # saveRDS(toRet, paste0("./generated_scTFinputs/LPS_",interpl,"_uncentered_",num_used_timepoints,"timepts.used.selected_",num_archetypes,"archetypes_nokmeans_pccomps",pc_comps,".rds"))
      # saveRDS(toRet, paste0("./generated_scTFinputs/LPS_",interpl,"_uncentered_",num_used_timepoints,"timepts.used_",num_archetypes,"archetypes_nokmeans_pccomps",pc_comps,".rds"))
      
    }
    
    #for original pca----
    data = macro_LPS@meta.data
    pca = data.frame(macro_LPS@reductions$pca@cell.embeddings)
    data= cbind(data, pca[ match(rownames(data), rownames(macro_LPS@reductions$pca)),])
    head(data)
    ggplot(data, aes(PC_1,PC_2, group=cell))+
      geom_path(aes(group = cell),size=0.1, alpha=0.25)+
      geom_point(aes(color = as.factor(timept_num)),size=1.5)+theme_classic(base_size = 16)
    
    #var explained-----
    t.data = macro_LPS@assays$RNA@data
    pca <- prcomp(t.data, scale = F, center =T)
    pca_scores = pca$x
    pca_loadings = pca$rotation
    sdev = pca$sdev
    var = unlist(sdev^2)
    pve = unlist(round(var/sum(var) * 100, 
                       digits = 6))
    pve
    plot(log10(pve[1:60]), pch = 19)
    sum((pve[1:50]))
    
    # pheatmap(cor(data[,10:59],data[,10:59]), cluster_rows = F, cluster_cols = F)
    # pheatmap(cor(data[,10:59], mean_pcscores[,c(4:53)]), cluster_rows = F, cluster_cols = F)
    # pheatmap(cor(data[,10:59], spline_pts_PC1.m[,c(3,5:53)]), cluster_rows = F, cluster_cols = F)
    
    #ground truth vs. random connection similarity between every pair of matrices-----
    library(lsa);library(proxy)
    collect_cosine=matrix(nrow = length(unique(data$cell)), ncol = length(unique(mean_pcscores$bin)))
    rownames(collect_cosine) = unique(data$cell)
    colnames(collect_cosine) = unique(mean_pcscores$bin)
    i=1
    for (cell in unique(data$cell)){
      A.subset = (data[cell== data$cell, 10:(9+num_comps)])
      j=1
      for (archetype in unique(mean_pcscores$bin)){
        print(archetype)
        B.subset = (mean_pcscores[archetype==mean_pcscores$bin, c(4:(3+num_comps))])
        
        cosine.sim = proxy::dist(A.subset, B.subset, pairwise=T, method="Euclidean")
        # cosine.sim = 1-proxy::simil(A.subset, B.subset, pairwise=T, method="cosine")
        collect_cosine[i, j] = sum(cosine.sim)
        # collect_cosine[i, j] = mean(cosine.sim)
        j=j+1
      }
      i=i+1
    }
    pheatmap(collect_cosine)
    new.df = data.frame(max.random = apply(collect_cosine, 1, max, na.rm=TRUE),
                        min.random = apply(collect_cosine, 1, min, na.rm=TRUE),
                        mean.random=rowMeans(collect_cosine))
    colnames(new.df) = paste0(colnames(new.df), "_pccomps", num_comps, "_", num_used_timepoints)
    metrics = cbind(metrics, new.df)
    
    
    #ground truth vs. scREALTIME (used 500 trajectories)-----
    collect_cosine=matrix(nrow = length(unique(data$cell)), ncol = length(unique(spline_pts_PC1.m$bin)))
    rownames(collect_cosine) = unique(data$cell)
    colnames(collect_cosine) = unique(spline_pts_PC1.m$bin)
    i=1
    for (cell in unique(data$cell)){
      print(cell)
      A.subset = (data[cell==data$cell, 10:(9+num_comps)])
      j=1
      for (archetype in unique(spline_pts_PC1.m$bin)){
        B.subset = (spline_pts_PC1.m[archetype==spline_pts_PC1.m$bin, c(3,5:(3+num_comps))])
        
        cosine.sim = proxy::dist(A.subset, B.subset, pairwise=T, method="Euclidean")
        # cosine.sim = 1-proxy::simil(A.subset, B.subset, pairwise=T, method="cosine")
        collect_cosine[i, j] = sum(cosine.sim)
        # collect_cosine[i, j] = mean(cosine.sim)
        j=j+1
      }
      i=i+1
    }
    pheatmap(collect_cosine)
    new.df = data.frame(max.scREAL = apply(collect_cosine, 1, max, na.rm=TRUE),
                        min.scREAL = apply(collect_cosine, 1, min, na.rm=TRUE),
                        mean.scREAL= rowMeans(collect_cosine))
    colnames(new.df) = paste0(colnames(new.df), "_pccomps", num_comps,"_", num_used_timepoints)
    metrics = cbind(metrics, new.df )
    # colnames(metrics)[26:31]=paste0(colnames(metrics[26:31]), ".selected")
    # colnames(metrics)[122:151]=paste0(colnames(metrics[122:151]), ".selected")
    metrics.all.m = melt(metrics)
    
    levels(metrics.all.m$variable) = gsub("random_min","min.random", levels(metrics.all.m$variable))
    metrics.all.m$timepts_used = gsub(".*_", "", metrics.all.m$variable)
    metrics.all.m$pccomps = gsub("pccomps","",sapply(strsplit(as.character(metrics.all.m$variable), "_"),"[[", 2) )
    ggplot(metrics.all.m[grepl("min",metrics.all.m$variable)&grepl("_5$|selected",metrics.all.m$variable),], aes(variable, value))+
      geom_jitter(width = 0.2,size=0.5 )+geom_boxplot(aes(color =timepts_used),outlier.shape = NA, alpha=0.05)+ 
      ylab("euclidean dist to ground truth")+
      theme_bw(base_size = 14)+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    
  }
  
}
metrics.all.m$pccomp.normalized = ifelse(metrics.all.m$pccomps==2, metrics.all.m$value/2,
                                         ifelse(metrics.all.m$pccomps==5, metrics.all.m$value/5,
                                                ifelse(metrics.all.m$pccomps==10, metrics.all.m$value/10,
                                                       ifelse(metrics.all.m$pccomps==20, metrics.all.m$value/20,
                                                              ifelse(metrics.all.m$pccomps==50, metrics.all.m$value/50,
                                                                     metrics.all.m$value)))))
metrics.all.m$value.normalized = ifelse(metrics.all.m$timepts_used==10, metrics.all.m$value/sqrt(2),
                                        ifelse(metrics.all.m$timepts_used==20, metrics.all.m$value/sqrt(4),
                                               ifelse(metrics.all.m$timepts_used==50, metrics.all.m$value/sqrt(10),
                                                      metrics.all.m$value)))
metrics.all.m$both.normalized = ifelse(metrics.all.m$timepts_used==10, metrics.all.m$pccomp.normalized/sqrt(2),
                                       ifelse(metrics.all.m$timepts_used==20, metrics.all.m$pccomp.normalized/sqrt(4),
                                              ifelse(metrics.all.m$timepts_used==50, metrics.all.m$pccomp.normalized/sqrt(10),
                                                     metrics.all.m$pccomp.normalized)))
ggplot(metrics.all.m[grepl("min",metrics.all.m$variable),], aes(as.factor(variable), log(both.normalized)))+
  geom_jitter(width = 0.2,size=0.5 )+geom_boxplot(aes(color =pccomps),outlier.shape = NA, alpha=0.05)+ 
  ylab("log(normalized euclidean dist to ground truth)")+ facet_grid(~timepts_used,scales="free")+
  theme_bw(base_size = 10)+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# metrics$metrics = rownames(metrics)
# write.table(metrics,"./generated_scTFinputs/distance_Euclidean_randomVSscREAL_paths.txt",quote=F,sep="\t",row.names = F)
# write.table(metrics,"./generated_scTFinputs/distance_Euclidean_randomVSscREAL_paths_addtxnnoise.txt",quote=F,sep="\t",row.names = F)
# write.table(metrics,"./generated_scTFinputs/distance_Euclidean_randomVSscREAL_paths_addtxnnoise10p.txt",quote=F,sep="\t",row.names = F)
# write.table(metrics,"./generated_scTFinputs/distance_Euclidean_randomVSscREAL_paths_otherPCcomps.txt",quote=F,sep="\t",row.names = F)
# write.table(metrics,"./generated_scTFinputs/distance_Euclidean_randomVSscREAL_paths_addtxnnoise10p_otherPCcomps.txt",quote=F,sep="\t",row.names = F)

# write.table(metrics,"./generated_scTFinputs/distance_CosineDist_randomVSscREAL_paths_otherPCcomps.txt",quote=F,sep="\t",row.names = F)
# write.table(metrics,"./generated_scTFinputs/distance_CosineDist_randomVSscREAL_paths_addtxnnoise10p_otherPCcomps.txt",quote=F,sep="\t",row.names = F)


# write.table(metrics,"./generated_scTFinputs/distance_Euclidean_randomVSscREAL_paths_otherPCcomps_fixed.txt",quote=F,sep="\t",row.names = F)
# write.table(metrics,"./generated_scTFinputs/distance_Euclidean_randomVSscREAL_paths_addtxnnoise10p_otherPCcomps_fixed.txt",quote=F,sep="\t",row.names = F)

# write.table(metrics,"./generated_scTFinputs/distance_CosineDist_randomVSscREAL_paths_otherPCcomps_fixed.txt",quote=F,sep="\t",row.names = F)
# write.table(metrics,"./generated_scTFinputs/distance_CosineDist_randomVSscREAL_paths_addtxnnoise10p_otherPCcomps_fixed.txt",quote=F,sep="\t",row.names = F)

metrics1 = read.delim("./generated_scTFinputs/distance_Euclidean_randomVSscREAL_paths_otherPCcomps.txt");metrics1$metrics="0p.txnnoise";
colnames(metrics1) = gsub("random_min","min.random", colnames(metrics1))
colnames(metrics1) = gsub("random_max","max.random", colnames(metrics1))
colnames(metrics1) = gsub("random_mean","mean.random", colnames(metrics1))
metrics2 = read.delim("./generated_scTFinputs/distance_Euclidean_randomVSscREAL_paths_addtxnnoise.txt");metrics2$metrics="2p.txnnoise"
metrics3 = read.delim("./generated_scTFinputs/distance_Euclidean_randomVSscREAL_paths_addtxnnoise10p_otherPCcomps.txt");metrics3$metrics="10p.txnnoise"

metrics.all = rbind(metrics1,metrics3) #metrics2, 
metrics.all.m = melt(metrics.all)
metrics.all.m$timepts_used = gsub(".*_", "", metrics.all.m$variable)
metrics.all.m$pccomps = gsub("pccomps","",sapply(strsplit(as.character(metrics.all.m$variable), "_"),"[[", 2) )
metrics.all.m$pccomp.normalized = ifelse(metrics.all.m$pccomps==2, metrics.all.m$value/sqrt(2),
                                         ifelse(metrics.all.m$pccomps==5, metrics.all.m$value/sqrt(5),
                                                ifelse(metrics.all.m$pccomps==10, metrics.all.m$value/sqrt(10),
                                                       ifelse(metrics.all.m$pccomps==20, metrics.all.m$value/sqrt(20),
                                                              ifelse(metrics.all.m$pccomps==50, metrics.all.m$value/sqrt(50),
                                                                     metrics.all.m$value)))))
metrics.all.m$value.normalized = ifelse(metrics.all.m$timepts_used==10, metrics.all.m$value/(10),
                                        ifelse(metrics.all.m$timepts_used==20, metrics.all.m$value/(20),
                                               ifelse(metrics.all.m$timepts_used==50, metrics.all.m$value/(50),
                                                      metrics.all.m$value/(5))))
metrics.all.m$both.normalized = ifelse(metrics.all.m$pccomps==2, metrics.all.m$value.normalized/sqrt(2),
                                       ifelse(metrics.all.m$pccomps==5, metrics.all.m$value.normalized/sqrt(5),
                                              ifelse(metrics.all.m$pccomps==10, metrics.all.m$value.normalized/sqrt(10),
                                                     ifelse(metrics.all.m$pccomps==20, metrics.all.m$value.normalized/sqrt(20),
                                                            ifelse(metrics.all.m$pccomps==50, metrics.all.m$value.normalized/sqrt(50),
                                                                   metrics.all.m$value.normalized)))))
metrics.all.m$timepts_used = factor(metrics.all.m$timepts_used, levels = c("5", "10","20","50","5.selected"))
ggplot(metrics.all.m[grepl("min",metrics.all.m$variable)&grepl("^50$",metrics.all.m$pccomps),], aes(variable, log(both.normalized)))+
  geom_jitter(width = 0.2,size=0.5)+geom_boxplot(aes(color =timepts_used),outlier.shape = NA, alpha=0.05)+ 
  ylab("log(euclidean dist to ground truth)")+
  facet_grid(~metrics, scales = "free")+
  theme_bw(base_size = 14)+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
metrics.all.m$pccomps = factor(metrics.all.m$pccomps, levels = c("2","5", "10","20","50"))
ggplot(metrics.all.m[grepl("min",metrics.all.m$variable)&grepl("5$|5.selected",metrics.all.m$timepts_used)&!grepl("^2$",metrics.all.m$pccomps),], aes(variable, log(value.normalized)))+
  geom_jitter(width = 0.2,size=0.5)+geom_boxplot(aes(color =pccomps),outlier.shape = NA, alpha=0.05)+ 
  ylab("log(euclidean dist to ground truth)")+
  facet_grid(~metrics+timepts_used, scales = "free")+
  theme_bw(base_size = 14)+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

#plot reconstruction----
stim = "LPS"
reduction_name = "pca"
interpl = "spline"
num_archetypes = 50
plots= list()
for (stim in c("LPS")){#c("LPS", "PIC","IFNb", "P3CSK","CpG","TNF")){
  print(stim)
  reconstructed_pc = readRDS(paste0("./generated_scTFinputs/LPS_",interpl,"_uncentered_",num_used_timepoints,"timepts.used.selected_",num_archetypes,"archetypes_nokmeans.rds"))
  reconstructed_pc.traj = (reconstructed_pc$reconstructed_trajectories)
  
  # reconstructed_pc.traj = (reconst$reconstructed_trajectories)
  reconstructed_pc.traj$mean = rowMeans(reconstructed_pc.traj[,1:(ncol(reconstructed_pc.traj)-2)]) 
  reconstructed_pc.traj$sdev = apply( (reconstructed_pc.traj[,1:(ncol(reconstructed_pc.traj)-3)]), 1, sd, na.rm=TRUE)
  p1=ggplot(reconstructed_pc.traj, aes(time, `AP1-n1-ktA0.48-kdeg30`)) + #ylim(0,0.3)+
    # geom_ribbon(aes(ymin = mean - sdev,
    #                 ymax = mean + sdev), alpha = 0.5)+
    # stat_summary(aes(time, `NFkB-n3-ktA0.96-kdeg30`),geom="line", color = "red",fun = "mean", linewidth=1, linetype="dashed")+
    geom_line(aes(group = as.factor(path), color = path), alpha = 0.25)+theme_bw(base_size = 8)+ggtitle(stim)
  p2=ggplot(reconstructed_pc.traj, aes(time, `IRF-n1-ktA0.48-kdeg30`)) + #ylim(0,0.1)+
    # geom_ribbon(aes(ymin = mean - sdev,
    #                 ymax = mean + sdev), alpha = 0.5)+
    # stat_summary(aes(time, `NFkB-n3-ktA0.96-kdeg30`),geom="line", color = "red",fun = "mean", linewidth=1, linetype="dashed")+
    geom_line(aes(group = as.factor(path), color = path), alpha = 0.25)+theme_bw(base_size = 8)+ggtitle(stim)
  p3=ggplot(reconstructed_pc.traj, aes(time, `NFkB-n1-ktA0.48-kdeg30`)) + #ylim(0,1)+
    # geom_ribbon(aes(ymin = mean - sdev,
    #                 ymax = mean + sdev), alpha = 0.5)+
    # stat_summary(aes(time, `NFkB-n3-ktA0.96-kdeg30`),geom="line", color = "red",fun = "mean", linewidth=1, linetype="dashed")+
    geom_line(aes(group = as.factor(path), color = path), alpha = 0.25)+theme_bw(base_size = 14)+ggtitle(stim)
  p4=ggplot(reconstructed_pc.traj, aes(time, `NFkB.IRF-n1-ktA0.48-kdeg30`)) + #ylim(0,0.1)+
    # geom_ribbon(aes(ymin = mean - sdev,
    #                 ymax = mean + sdev), alpha = 0.5)+
    # stat_summary(aes(time, `NFkB-n3-ktA0.96-kdeg30`),geom="line", color = "red",fun = "mean", linewidth=1, linetype="dashed")+
    geom_line(aes(group = as.factor(path), color = path), alpha = 0.25)+theme_bw(base_size = 8)+ggtitle(stim)
  
  p5=ggplot(reconstructed_pc.traj, aes(time, `NFkB.p38-n1-ktA0.48-kdeg30`)) + #ylim(0,0.1)+
    # geom_ribbon(aes(ymin = mean - sdev,
    #                 ymax = mean + sdev), alpha = 0.5)+
    # stat_summary(aes(time, `NFkB-n3-ktA0.96-kdeg30`),geom="line", color = "red",fun = "mean", linewidth=1, linetype="dashed")+
    geom_line(aes(group = as.factor(path), color = path), alpha = 0.25)+theme_bw(base_size = 8)+ggtitle(stim)
  p1/p2/p3/p4/p5
  
  ggplot(reconstructed_pc.traj, aes(time, `IRF-n1-ktA0.48-kdeg30`)) + #ylim(0,0.5)+
    # geom_ribbon(aes(ymin = mean - sdev,
    #                 ymax = mean + sdev), alpha = 0.5)+
    # stat_summary(aes(time, `NFkB-n3-ktA0.96-kdeg30`),geom="line", color = "red",fun = "mean", linewidth=1, linetype="dashed")+
    geom_line(aes(group = as.factor(path)),color = "forestgreen", alpha = 0.3)+theme_bw(base_size = 14)+ggtitle(stim)
  
  
  plots[[stim]]=NA
  ggplot((reconstructed_pc$reconstructed_trajectories), aes(time, `NFkB-n3-ktA0.96-kdeg30`)) + #ylim(0,10)+
    stat_summary(aes(time, `NFkB-n3-ktA0.96-kdeg30`),geom="line", color = "red",fun = "mean", linewidth=1, linetype="dashed")+
    geom_line(aes(group = as.factor(path), color = path), alpha = 0.25)+theme_bw(base_size = 8)+ggtitle(stim)
  
}



#get trajectory features for ground truth----------------------------------------
data = readRDS("./generated_scTFinputs/output_test_5GRS_all40paramsets_LPS.rds")

mat_allstims <- list()
num_paths_sum = 0
labels = c()
for (stim in c('LPS')){ #},'TNF1','PIC1','P3K1','CpG1')){
  print(stim)
  simulated.data = readRDS(paste0("./generated_scTFinputs/output_test_5GRS_all40paramsets_",stim,".rds"))
  
  num_paths = length(unique(simulated.data$cell))
  num_timepts = length(unique(simulated.data$time))
  
  mat_allstims <- rbind(mat_allstims, simulated.data)
  
  num_paths_sum = num_paths_sum + num_paths
  labels = c(labels, rep(stim, num_paths))
}
mat.numbers = mat_allstims[,!grepl("time|cell|label", colnames(mat_allstims))]
mat.meta = mat_allstims[,grepl("time|cell|label", colnames(mat_allstims))]
if (1){ #rescale 0-1 or not----
  mat.numbers = apply(mat.numbers, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X))) #rescale each gene column 0-1 over all stims
  # mat.numbers = cbind(mat.numbers, mat.meta)
} 

###############################################################################
# Figure 3: plot macrophage measured data ----
###############################################################################

# plot violins, group by timepoint----
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
# Figure 3: impute all trajectories ----
# Note: from F://BACKUP.../Projects.../trajectory_method/tensor_trajectory.R

###############################################################################

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
                           varFilter = T, exp_prob = 1) 
  
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

# for M2(IL4) macrophages------
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
# plot scaled across all stimuli, comment/uncomment for M0, M1, M2----
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

################################################################################
# plot scaled across all stimuli and mac types for comparing polarization states-------
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


###############################################################################
# plot trajectories----
mat.numbers2 = cbind(mat.numbers, mat.meta)
mat.numbers2$path_stim = paste0(mat.numbers2$path,"_", mat.numbers2$stimulus)
ggplot(mat.numbers2[grepl("",mat.numbers2$stimulus),], aes(time,Cxcl10, group = path_stim)) + ylim(0,1)+
  geom_vline(xintercept = c(0,0.25,1,3,8), linetype="dotted")+ #theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  # geom_vline(xintercept = c(0,0.5,1,3,5,8), linetype="dotted")+
  geom_line(aes(group = as.factor(path_stim), color = stimulus), alpha = 0.05)+
  theme_bw(base_size = 14)+theme(legend.position = "none")

mat.numbers2$path_stim = paste0(mat.numbers2$path,"_", mat.numbers2$stimulus,"_",mat.numbers2$type)
ggplot(mat.numbers2[grepl("LPS",mat.numbers2$stimulus)&grepl("",mat.numbers2$type),], aes(time,Nfkbia, group = path_stim)) + ylim(0,1)+
  # geom_vline(xintercept = c(0,0.25,1,3,8), linetype="dotted")+ #theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
  geom_vline(xintercept = c(0,0.5,1,3,5,8), linetype="dotted")+
  geom_line(aes(group = as.factor(path_stim), color = type), alpha = 0.05)+
  theme_bw(base_size = 14)+theme(legend.position = "none")
# single gene------
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


# multiple genes-----
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
# Figure 4: calculate trajectory features----
###############################################################################

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
# Figure 4: calculate mutual information on trajectory features----
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

###############################################################################
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

# make correlation violin/bar plots------
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




###############################################################################
# Figure 6: polarized macrophages
###############################################################################
####################################################### plot the trajectories-----
# plot scaled across all stimuli, comment/uncomment for M0, M1, M2----
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
#plot PCA loadings-------
feature = "integral"#"speed1hr" #"peak_amp_lfc"
loadings = read.delim(paste0("./trajectory/matrix_trajectory_features_M0M1M2_",feature,"_k20_prcomp_loadings.txt"))
scores = read.delim(paste0("./trajectory/matrix_trajectory_features_M0M1M2_",feature,"_k20_prcomp_scores.txt"))
clusters = readxl::read_excel("F://scRNAseq_macro/SuppTables/TableS4_gene_regulatory_strategies_allgenes.xlsx")
loadings$GRS = clusters$clusters[match(loadings$Loading, clusters$gene)]
tmp = loadings[,c(1:3,496)]
colors_list = c("darkorange", "gold","darkgreen", "lightblue","darkred")
colors_list = c("darkorange", "forestgreen","red", "blue1","darkred")
ggplot(loadings, aes(PC1, PC2))+geom_point(aes(fill = GRS), pch=21, size=3, alpha=0.8)+
  scale_fill_manual(values = colors_list)+
  theme_bw(base_size = 16)+ggtitle(feature)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)

#collect dim best start with 1D top 20 (done on Precision)---- 

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



###############################################################################
# Figure 7: Response Dynamics allows better polarization state identification----
###############################################################################

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

