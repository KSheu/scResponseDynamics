# single-cell response trajectory imputation algorithm
# auxilliary functions for plotting and data extraction
# contact: katherinesheu[at]ucla[dot]edu


#' scREALTIME
#' Input a Seurat object of single cell data and returns list object containing single cell trajectories
#' @param macro Seurat object of single cell data with multiple measured timepoints
#' @param metadata Metadata frame of the measured single cells, obtain with getMetaData
#' @param num_archetypes = 20 Number of cell archetypes expected per timepoint (eg. 10, 20)
#' @param timepoints Timepoints to use from the measured data
#' @param num_trajectories = 1000 Number of simulated trajectories expected
#' @param num_sim_pts = 100 Number of simulated timepoints to output
#' @param reduction = "pca" Dimensionality reduction to use (eg. "pca","ica","nmf")
#' @param stimulus Specify one stimulus to impute if desired
#' @param consensus_measure Method to define cell archetypes ("mean", "median")
#' @param interpolant Method to interpolate timepoints ("spline","linear")
#' @param data = "RNA" normalization from Seurat object to use (eg. "SCT", "ISNorm")
#' @param prob_method = "distance" Options: 'distance', 'density', 'hybrid'. Linkage probability based on a 'distance' metric, 'density' of cells in archetypes, or 'hybrid'
#' @param distance_metric = "euclidean" Distance metric to use to identify cell archetype links over timepoints. This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param varFilter = T Filters out zero variance genes prior to kmeans clustering.
#' @param exp_prob = 1 Raises the transition probability matrix to the power of exp_prob. Higher values make weak links weaker.
#' @return trajectoryObject

scREALTIME = function(macro, metadata, num_archetypes=20, timepoints, num_trajectories = 1000, num_sim_pts = 100,
                      reduction = "pca", stimulus, consensus_measure = "median", interpolant="spline", data = "RNA", prob_method = 'distance', distance_metric = 'euclidean', varFilter = T, exp_prob = 1){


  require(NbClust); require(matrixStats); require(Rfast);require(Seurat);require(devtools)
  require(RColorBrewer); require(ggplot2); require(gridExtra);


  # retrieve desired data subset------------------------------------------------
  cells_by_timept <- list()
  for(i in timepoints){
    index = paste("time_",i, "hr", sep = "")
    query = paste(i, "hr", sep='')
    #cells_by_timept[[index]] <- rownames(metadata)[metadata$timept == query]
    cells_by_timept[[index]] <- rownames(metadata)[metadata$timept_num == i]
  }

  if(data == 'RNA'){
    RNA <- as.data.frame(macro@assays$RNA@data)
  }
  else if(data == 'ISNorm'){
    RNA <- as.data.frame(macro@assays$ISnorm@data)
  }

  RNA <- t(RNA)

  # if(toupper(reduction) == 'pca'){
  #   pca = prcomp(RNA, center = F, scale = F, rank. = 50)
  #   pcscores = pca$x
  #   loadings = pca$rotation
  #   # pcscores = macro[['pca']]@cell.embeddings
  #   #print(pcscores)
  # }else if(toupper(reduction) == 'NMF'){
  #   pcscores = macro[['NMF']]@cell.embeddings
  # }else if(toupper(reduction) == 'ICA'){
  #   pcscores = macro[['ica']]@cell.embeddings
  # }

  RNA <- RNA[rownames(RNA) %in% rownames(metadata),]

  clusterings <- list()
  zero_var_inds <- list()

  # kmeans clustering ----------------------------------------------------------
  for(i in timepoints){
    index = paste("time_", i, "hr", sep = "")
    data <- RNA[rownames(RNA) %in% cells_by_timept[[index]],]
    zero_var_inds[[index]] <- resample::colVars(data) == 0
    print(zero_var_inds)
    print(tail(data[,! zero_var_inds[[index]]]))
    print(dim(data[,! zero_var_inds[[index]]]))

    if(varFilter){
      clusterings[[index]] <- kmeans(data[,!zero_var_inds[[index]]], centers = num_archetypes, iter.max = 50)
    }else{
      clusterings[[index]] <- kmeans(data, centers = num_archetypes, iter.max = 50)
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

  # Plotting
  # clustering_plots <- list()
  # for(i in timepoints){
  #   index = paste("time_", i, "hr", sep = "")
  #   clustering_plots[[index]] <- fviz_cluster(clusterings[[index]], data = RNA[rownames(RNA) %in% cells_by_timept[[index]],!zero_var_inds[[index]]],
  #                                             geom = "point",
  #                                             ellipse.type = "convex",
  #                                             ggtheme = theme_bw(), main = FALSE
  #   )
  #
  # }
  # ggpubr::ggarrange(plotlist = clustering_plots, labels = paste(timepoints, 'hr'), hjust = -1.75, vjust = 1.75)
  #
  # Here we will use our own PCA
  #pcscores <- macro[[reduction]]@cell.embeddings
  #df_for_pca = data.frame(RNA, metadata$timept)
  #colnames(df_for_pca)
  #View(df_for_pca)

  # retrieve Seurat object PCA (or other dim reduction)--------------------------
  if(toupper(reduction) == 'PCA'){
    # pca = prcomp(RNA, center = F, scale = F, rank. = 50)
    # pcscores = pca$x
    # loadings = pca$rotation
    pcscores = macro[['pca']]@cell.embeddings
    #print(pcscores)
  }else if(toupper(reduction) == 'NMF'){
    pcscores = macro[['NMF']]@cell.embeddings
  }else if(toupper(reduction) == 'ICA'){
    pcscores = macro[['ica']]@cell.embeddings
  }


  cluster_densities <- as.data.frame(prop.table(as.matrix(cluster_counts), 2))
  cluster_densities <- cluster_densities^(exp_prob)

  # identify cell archetypes in PC space ----------------------------------------
  cell_cluster_df <- matrix(nrow = 0, ncol = 2)
  for(i in timepoints){
    index = paste("time_", i, "hr", sep = "")
    cell_cluster_df <- c(cell_cluster_df, clusterings[[index]]$cluster)
  }
  cell_cluster_df <- as.data.frame(cell_cluster_df)


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
  #View(mean_pcscores)
  #print(str(mean_pcscores))

  # Generating walk probabilities matrix -----------------------------------------
  distances = as.matrix(dist(mean_pcscores[,3:ncol(mean_pcscores)],method = distance_metric))
  rownames(distances) = mean_pcscores$Group.1
  colnames(distances) = mean_pcscores$Group.1

  # Normalize so larger distances have lower probability

  walks <- list()
  walk_probs <- matrix(nrow = num_trajectories, ncol = 2)
  walk_probs[,1] = c(1:num_trajectories)

  if(prob_method == 'distance'){
    print("Generating random walks based on distance between clusters")
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


  # Spline fitting across the random walks------------------------------------------

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

  # Now, we get number of unique paths * number of trajectories * 50 PCs
  num_splines =dim(walk_frequencies)[1] * 50
  print(paste("Number of splines: ", dim(walk_frequencies)[1] * 50))
  dup_rows = duplicated(spline_pts[[1]])
  walk_probs = walk_probs[!dup_rows, ]

  for(i in 1:length(spline_pts)){
    spline_pts[[i]] <- spline_pts[[i]][!dup_rows,]

    #spline_pts[[i]] <- unique(spline_pts[[i]])
  }

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

  simulated = matrix(ncol = (ncol(mean_pcscores)-3), nrow = num_sim_pts*dim(walk_frequencies)[1])
  simulation_mat_time = rep(sim_times, dim(walk_frequencies)[1])

  for(i in 1:(ncol(mean_pcscores)-3)){ # Loop over PCs
    preds = c()
    for(walk in 1:nrow(spline_pts[[i]])){
      if(interpolant == 'spline'){
        spline <- smooth.spline(timepoints, spline_pts[[i]][walk,], cv = T, all.knots = T) # TIMEPOINTS
        preds <- c(preds, predict(spline, sim_times)$y)
      } else if(interpolant == 'linear'){
        fun = approxfun(timepoints, spline_pts[[i]][walk,], rule = 2)
        preds = c(preds, fun(sim_times))
      } else if(interpolant == 'loess'){
        loe = loess(spline_pts[[i]][walk,] ~ timepoints, span = 0.50, se = F)
        preds = c(preds, predict(loe, sim_times))
      }else {
        print('Unspecified interpolant')
      }
    }
    simulated[,i] = preds
  }

  # Cast back to gene expression space ------------------------------------------------
  loadings <- macro[[reduction]]@feature.loadings

  if(toupper(reduction) == 'NMF'){
    loadings = macro[['NMF']]@feature.loadings
  }
  if(toupper(reduction) == 'ICA'){
    loadings = macro[['ica']]@feature.loadings
  }
  #dim(loadings)
  reconstructed_pc <- t(loadings %*% t(simulated))
  reconstructed_pc <- as.data.frame(reconstructed_pc)
  reconstructed_pc$time <- simulation_mat_time
  reconstructed_pc$path <- rep(1:dim(walk_frequencies)[1], each = num_sim_pts)

  # return result --------------------------------------------------------------------
  toRet = list()
  toRet[['reconstructed_trajectories']] = reconstructed_pc
  toRet[['cluster_densities']] = cluster_densities
  toRet[['metadata']] = metadata
  toRet[['number_unique_trajectories']] = num_unique_traj
  toRet[['number_of_splines']] = num_splines
  toRet[['probability_of_walks']] = walk_probs

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

  return(toRet)
}

#' Input a Seurat object of single cell data and returns list object containing single cell trajectories
#' @param macro Seurat object of single cell data with multiple measured timepoints
#' @param stimulus Stimulus name to label/subset the dataset
#' @param timepoints Vector of numeric timepoints in the dataset

getMetaData = function(macro, stimulus, timepoints){
  require('Seurat')
  require('factoextra')
  require('matrixStats')
  require('cowplot')
  require('ggplot2')
  require('RColorBrewer')
  require('NbClust')
  require('devtools')
  #source('ks_scripts/ksheu.library1/R/NNMF_from_file.R')
  require('NMF')
  require('gridExtra')
  require('Rfast')
  pca <- macro[['pca']]
  #PCAplot(macro)
  metadata <- as.data.frame(macro@meta.data)

  order = paste(timepoints,'hr', sep = '')
  #View(head(as.data.frame(macro@assays$RNA@data)))
  # Write file for NNMF
  #nnmf=nmf(as.data.frame(macro@assays$RNA@data), rank = 5)

  plotname = paste('PCA Plot of', stimulus, 'stimulated data')
  # Confine our metadata to stimulus and Unstim 0 hr
  metadata <- metadata[metadata$stimulus == stimulus | metadata$timept == '0hr'| metadata$timept == '0.0hr',]
  # print(PCAPlot(macro, cells = as.vector(rownames(metadata)), group.by = 'timept', order = order, shuffle = TRUE, label = TRUE, label.box = TRUE) + ggtitle(plotname))
  # See clustering based on time point, but clusters vary enough, that perhaps we can subdivide them into archetypes
  # For flexibility with number of clusters, we can use k-means clustering where k = number of archetypes

  metadata$timept_num = as.numeric(sapply(strsplit(metadata$timept,"h"), `[`, 1))

  metadata = metadata[trunc(metadata$timept_num, 5) %in% trunc(timepoints, 5),]

  return(metadata)
}
