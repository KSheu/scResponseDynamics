# single-cell response trajectory imputation algorithm
# auxilliary functions for plotting and data extraction
# contact: katherinesheu[at]ucla[dot]edu


#' getTrajectory4
#' Input a Seurat object of single cell data and returns list object containing single cell trajectories
#' @param macro Seurat object of single cell data with multiple measured timepoints
#' @param metadata Metadata frame of the measured single cells, obtain with getMetaData
#' @param num_archetypes Number of cell archetypes expected per timepoint (eg. 10, 20)
#' @param timepoints Timepoints to use from the measured data
#' @param num_trajectories Number of simulated trajectories expected
#' @param num_sim_pts Number of simulated timepoints to output
#' @param reduction Dimensionality reduction to use (eg. "pca","ica","nmf")
#' @param stimulus Specify one stimulus to impute if desired
#' @param consensus_measure Method to define cell archetypes ("mean", "median")
#' @param interpolant Method to interpolate timepoints ("spline","linear")
#' @param data = "RNA" normalization from Seurat object to use (eg. "SCT", "ISNorm")
#' @param prob_method = 'distance'
#' @param distance_metric = 'euclidean'
#' @param varFilter = T
#' @param exp_prob = 1
#'  

getTrajectory4 = function(macro, metadata, num_archetypes, timepoints, num_trajectories, num_sim_pts, reduction, stimulus, consensus_measure, interpolant, data = "RNA", prob_method = 'distance', distance_metric = 'euclidean', varFilter = T, exp_prob = 1){
  
  cells_by_timept <- list()
  for(i in timepoints){
    index = paste("time_",i, "hr", sep = "")
    query = paste(i, "hr", sep='')
    #cells_by_timept[[index]] <- rownames(metadata)[metadata$timept == query]
    cells_by_timept[[index]] <- rownames(metadata)[metadata$timept_num == i]
  }
  # View(cells_by_timept)
  
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

  
  #plot for sanity
  require('factoextra')
  require('ggfortify')
  
  #fviz_eig(pca)
  #autoplot(pca, data = df_for_pca, loadings = FALSE, colour = 'metadata.timept')
  
  cluster_densities <- as.data.frame(prop.table(as.matrix(cluster_counts), 2))
  cluster_densities <- cluster_densities^(exp_prob)
  
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
  
  # Generating walk probabilities matrix
  distances = as.matrix(dist(mean_pcscores[,3:ncol(mean_pcscores)]))
  rownames(distances) = mean_pcscores$Group.1
  colnames(distances) = mean_pcscores$Group.1
  
  # Need to convert this to a transitions probs matrix at each step
  
  # Need to normalize so larger distances have lower probability
  
  
  
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
  
  # Now let's cast back to gene expression space
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
  # 
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
  print(PCAPlot(macro, cells = as.vector(rownames(metadata)), group.by = 'timept', order = order, shuffle = TRUE, label = TRUE, label.box = TRUE) + ggtitle(plotname))
  # See clustering based on time point, but clusters vary enough, that perhaps we can subdivide them into archetypes
  # For flexibility with number of clusters, we can use k-means clustering where k = number of archetypes
  
  metadata$timept_num = as.numeric(sapply(strsplit(metadata$timept,"h"), `[`, 1))
  
  metadata = metadata[trunc(metadata$timept_num, 5) %in% trunc(timepoints, 5),]
  
  return(metadata)
}

getElbowPlots = function(macro, metadata, stimulus, timepoints, k.max = 15, varFilter = T){
  require('matrixStats')
  require('purrr')
  require('ggplot2')
  plotlist = list()
  zero_var_inds <- list()
  RNA <- as.data.frame(macro@assays$RNA@data)
  RNA <- t(RNA)
  RNA <- RNA[rownames(RNA) %in% rownames(metadata),]
  
  cells_by_timept <- list()
  for(i in timepoints){
    index = paste("time_",i, "hr", sep = "")
    query = paste(i, "hr", sep='')
    cells_by_timept[[index]] <- rownames(metadata)[metadata$timept == query]
  }
  
  
  for(t in timepoints){
    index = paste("time_",t,"hr", sep = "")
    data <- RNA[rownames(RNA) %in% cells_by_timept[[index]],]
    zero_var_inds[[index]] <- colVars(data) == 0
    
    
    
    if(varFilter){
      # Use map_dbl to run many models with varying value of k (centers)
      tot_withinss <- map_dbl(1:k.max,  function(k){
        model <- kmeans(x = data[,!zero_var_inds[[index]]], centers = k, iter.max = 25)
        model$tot.withinss
      })
    }else{
      # Use map_dbl to run many models with varying value of k (centers)
      tot_withinss <- map_dbl(1:k.max,  function(k){
        model <- kmeans(x = scaled_data, centers = k, iter.max = 25)
        model$tot.withinss
      })
    }
    
    
    # Generate a data frame containing both k and tot_withinss
    elbow_df <- data.frame(
      k = 1:k.max,
      tot_withinss = tot_withinss
    )
    
    
    # Plot the elbow plot
    plotlist[[index]] = ggplot(elbow_df, aes(x = k, y = tot_withinss)) +
      geom_line() + geom_point()+
      scale_x_continuous(breaks = 1:k.max)
  }
  
  ggpubr::ggarrange(plotlist = plotlist, labels = paste(timepoints, 'hr'), font.label = list(size = 7, color = "blue"))
  
  
}

getSummaryStats = function(macro, metadata, reconstructed_data, gene, timepoints, scale, omitzero, timewindow=1){
  require('matrixStats')
  toRet = list()
  
  gene_ind = grep(pattern = gene, x = colnames(reconstructed_data[['reconstructed_trajectories']]))
  candidate_genes = colnames(reconstructed_data[['reconstructed_trajectories']])[gene_ind]
  #print(gene_ind)
  #print(colnames(reconstructed_data[['reconstructed_trajectories']])[gene_ind])
  if(length(gene_ind)==0){
    print('Gene not found in simulated data')
    return('Gene not found')
  }
  else{
    if(gene %in% candidate_genes){
      print(which(candidate_genes == gene))
      gene_ind = gene_ind[which(candidate_genes == gene)]
    }
    else{
      gene_ind = gene_ind[1] # uses first instance of desired gene
    }
  }
  test_gene1 = reconstructed_data[['reconstructed_trajectories']][,c(gene_ind,ncol(reconstructed_data$reconstructed_trajectories)-2, ncol(reconstructed_data$reconstructed_trajectories)-1)]
  #View(test_gene1)
  colnames(test_gene1)[1] = 'exp'
  if(scale){
    test_gene1$exp = (test_gene1$exp - min(test_gene1$exp))/(max(test_gene1$exp)-min(test_gene1$exp))
  }
  
  truedata = as.data.frame(macro@assays$RNA@data)
  gene_ind_true = grep(pattern = gene, x = rownames(truedata))
  if(length(gene_ind_true)==0){
    print('Gene not found in true data')
    return('Gene not found')
  }
  else{
    gene_ind_true = gene_ind_true[1] # uses first instance of desired gene
  }
  # print('these are cols')
  # print(colnames(truedata))
  # 
  # print('these are rows')
  # print(rownames(truedata))
  truedata = truedata[gene_ind_true,]
  if(omitzero){
    #truedata = truedata[,truedata != 0]
    truedata = truedata[,-(which(colSums(truedata)==0))]
  }
  #print(dim(truedata))
  if(is.null(dim(truedata)) || 0 %in% dim(truedata)){
    print(paste("the following gene has been omitted due to insufficient data after filtering:", gene))
    return(NULL)
  }
  #print(paste('dimension of filtered data', dim(truedata)))
  #print(truedata)
  if(scale){
    truedata[1,] = (truedata[1,]-min(truedata[1,]))/(max(truedata[1,])- min(truedata[1,]))
  }
  #View(truedata)
  #View(test_gene1)
  sim_means = c()
  sim_var = c()
  data_means = c()
  data_var = c()
  boxplots = list()
  vlnplots = list()
  ks = list()
  
  for(time in timepoints){
    timeptstr = paste(time, 'hr', sep = "")
    print(paste('Checking time', timeptstr))
    
    min_time = test_gene1$time[which.min(abs(time-(timewindow/2)-test_gene1$time))]
    max_time = test_gene1$time[which.min(abs(time+(timewindow/2)-test_gene1$time))]
    times = unique(test_gene1$time[test_gene1$time > min_time & test_gene1$time < max_time])
    #print(times)
    
    relevant_cells = rownames(metadata)[metadata$timept_num == time] # Get all cells that were measured at given time point 
    #print(length(relevant_cells))
    #print(relevant_cells)
    relevant_cells =  colnames(truedata) %in% relevant_cells
    #print(length(relevant_cells))
    #print(length(colnames(truedata)))
    #print('here is the data')
    #View(truedata[1,relevant_cells])
    truecells = cbind(t(truedata[1,relevant_cells]), 'True_Data')
    simcells = cbind(test_gene1$exp[test_gene1$time %in% times], 'Simulated_Data')
    plottabledf = as.data.frame(rbind(truecells, simcells))
    colnames(plottabledf) = c('value', 'group')
    plottabledf$group = as.factor(plottabledf$group)
    
    #View(plottabledf)
    #print(str(plottabledf))
    if(length(truecells[,1] !=0) & length(simcells[,1]) !=0){
      ks[[timeptstr]]=ks.test(unlist(as.numeric(truecells[,1])), unlist(as.numeric(simcells[,1])))
    }else{
      ks[[timeptstr]] == 'KS test could not run due at this time point to insufficient data'
    }

    #View(truecells)
    #View(simcells)
    
    boxplots[[timeptstr]] = ggplot(plottabledf, aes(x = group, y = value, fill = group)) + geom_boxplot() + ggtitle(paste('Comparison of simulated \n and true distributions of', gene, 'at', timeptstr)) + ylab(ifelse(scale, "scaled value", "unscaled value")) + geom_jitter(alpha=0.2)
    
    vlnplots[[timeptstr]] = ggplot(plottabledf, aes(x = group, y = value, fill = group)) + geom_violin() + ggtitle(paste('Comparison of simulated \n and true distributions of', gene, 'at', timeptstr)) + ylab(ifelse(scale, "scaled value", "unscaled value"))
    
    #View(as.matrix(truedata[1, relevant_cells %in% colnames(truedata)]))
    #datamean_cur = as.numeric(rowMeans(as.matrix(truedata[1, relevant_cells]))) # Get mean of true data for all relevant cells
    datamean_cur = as.numeric(mean(plottabledf$value[plottabledf$group == 'True_Data']))
    datavar_cur = as.numeric(var(plottabledf$value[plottabledf$group == 'True_Data']))
    #datavar_cur = as.numeric(rowVars(as.matrix(truedata[1, relevant_cells]))) # Get variance of true data for all relevant cells
    data_means = c(data_means, datamean_cur)
    data_var = c(data_var, datavar_cur)
    #print(data_means)
    
    #sim_means = c(sim_means, mean(test_gene1$exp[test_gene1$time == time]))
    sim_means = c(sim_means, as.numeric(mean(plottabledf$value[plottabledf$group == 'Simulated_Data'])))
    sim_var = c(sim_var, as.numeric(var(plottabledf$value[plottabledf$group == 'Simulated_Data'])))
    #print(sim_means)
    #sim_var = c(sim_var, var(test_gene1$exp[test_gene1$time == time]))
  }
  
  ret = cbind(paste(timepoints, 'hrs', sep = ''), data_means)
  ret = cbind(ret, data_var)
  ret = cbind(ret, sim_means)
  ret = cbind(ret, sim_var)
  colnames(ret) = c('Timepoint', paste('True mean', gene, sep = '-'), paste('True variance', gene, sep = '-'), paste('Simulated mean', gene, sep = '-'), paste('Simulated variances', gene, sep = '-'))
  
  ret = as.data.frame(ret)
  ret[,2:5] <- sapply(ret[,2:5],as.numeric)
  
  #View(truedata)
  #View(test_gene1)
  
  toRet[['stats']] = ret
  toRet[['boxplots']] = boxplots
  toRet[['gene']] = gene
  
  currenttimept = 1
  temp = list()
  for(p in boxplots){
    t = paste(timepoints[currenttimept], 'hr')
    temp[[t]] = p + ggtitle(t)
    currenttimept = currenttimept+ 1
  }
  
  toRet[['all_boxplots']] = ggpubr::annotate_figure(ggpubr::ggarrange(plotlist = temp, common.legend = TRUE, legend = "bottom"), top = paste("Comparison of true and simulated data for", gene))
  
  currenttimept = 1
  temp = list()
  for(p in vlnplots){
    t = paste(timepoints[currenttimept], 'hr')
    temp[[t]] = p + ggtitle(t)
    currenttimept = currenttimept+ 1
  }
  
  toRet[['all_vlnplots']] = ggpubr::annotate_figure(ggpubr::ggarrange(plotlist = temp, common.legend = TRUE, legend = "bottom"), top = paste("Comparison of true and simulated data for", gene))
  
  #View(plottabledf)
  toRet[['KStest']]=ks
  #toRet[['KStest']] = ks.test(unlist(test_gene1$exp), unlist(truedata[1,]))
  
  
  return(toRet)
}

generatePlots = function(reconstructed_data, gene, stimulus, scale=FALSE){
  
  
  num_sim_pts = 100
  toRet = list()
  gene_ind = grep(pattern = gene, x = colnames(reconstructed_data[['reconstructed_trajectories']]))
  print(gene_ind)
  candidate_genes = colnames(reconstructed_data[['reconstructed_trajectories']])[gene_ind]
  print(candidate_genes)
  if(length(gene_ind)==0){
    print('Gene not found in data')
    return('Gene not found')
  }
  else{
    if(gene %in% candidate_genes){
      print(which(candidate_genes == gene))
      gene_ind = gene_ind[which(candidate_genes == gene)]
    }
    else{
      gene_ind = gene_ind[1] # uses first instance of desired gene
    }
  }
  
  last2 = c(ncol(reconstructed_data[['reconstructed_trajectories']]), ncol(reconstructed_data[['reconstructed_trajectories']])-1)
  test_gene1 = reconstructed_data[['reconstructed_trajectories']][,c(gene_ind,last2)]
  test_gene1_name = colnames(test_gene1)[1]
  colnames(test_gene1)[1] <- 'exp'
  toRet[['unscaled_gene_exp']] = test_gene1$exp
  if(scale == 'max'){
    maxval = max(test_gene1$exp)
    test_gene1$exp = test_gene1$exp/maxval
  }
  else if(scale == 'zero-one'){
    maxval = max(test_gene1$exp)
    minval = min(test_gene1$exp)
    test_gene1$exp = (test_gene1$exp - minval)/(maxval-minval)
  }else if(scale == 'negatives2zero'){
    test_gene1$exp[test_gene1$exp<0] = 0
  }
  
  if(scale != 'none'){
    ylabel = 'scaled_exp'
  }
  else{
    ylabel = 'unscaled_exp'
  }
  
  # testing probability filter
  # high_prob_paths = reconstructed_data$probability_of_walks[reconst$probability_of_walks[,2] > 3e-5,1]
  # test_gene1 = test_gene1[test_gene1$path %in% high_prob_paths,]
  # 
  
  # plot_points = ggplot(data = test_gene1, aes(x = time, y = exp)) + geom_point(aes(color = factor(path))) + 
  #   xlab('time') + ylab(ylabel) + labs(title = paste(stimulus, '_', test_gene1_name)) + guides(color = FALSE) +
  #   scale_color_discrete()  + geom_smooth(method = 'gam')+theme_bw()
  
  plot_lines = ggplot(data = test_gene1, aes(x = time, y = exp, group = path)) + 
    geom_line(aes(color = factor(path)), alpha =0.2, size =0.1) +  xlab('time') + ylab(ylabel) + 
    labs(title = paste(stimulus, '_', test_gene1_name)) + guides(color = FALSE) +  
    scale_color_discrete()+theme_bw()
  
  # plot_few = ggplot(data = test_gene1[c(1:num_sim_pts*15),], aes(x = time, y = exp, group = path)) + 
  #   geom_line(aes(color = factor(path))) +  xlab('time') + ylab(ylabel) + 
  #   labs(title = paste(stimulus,'_', test_gene1_name)) + guides(color = FALSE) +  
  #   scale_color_discrete()+theme_bw()
  # 
  # plot_density = ggplot(data = test_gene1, aes(x = time, y = exp)) + 
  #   xlab('time') + ylab(ylabel) + 
  #   labs(title = paste(stimulus, '_', test_gene1_name)) + guides(color = FALSE) +  
  #   theme_classic()+
  #   # stat_density_2d(aes(fill = ..level..), geom="polygon")+
  #   geom_density_2d_filled(alpha = 1) +
  #   scale_fill_brewer(palette = "Greys")+ theme(legend.position = "none")+
  #   geom_line(aes(group = path), alpha = 0.2, position = position_dodge(.2)) 
  # 
  
  # plot_lines_sd = ggplot(data = test_gene1, aes(x = time, y = exp, group = path)) +
  #   geom_line(aes(x = time, y = exp, group = path), alpha =0.2, size =0.1)+
  #   xlab('time') + ylab(ylabel) +
  #   stat_summary(fun.y=median, fun.ymin=function(x) median(x) - sd(x)*3,
  #                fun.ymax=function(x) median(x) + sd(x)*3, geom="ribbon", fill = "green", alpha =0.2) +
  #   stat_summary(fun.y=median, fun.ymin=function(x) median(x) - sd(x)*2,
  #                fun.ymax=function(x) median(x) + sd(x)*2, geom="ribbon", fill = "green", alpha =0.2) +
  #   stat_summary(fun.y=median, fun.ymin=function(x) median(x) - sd(x),
  #                fun.ymax=function(x) median(x) + sd(x), geom="ribbon", fill = "green", alpha =0.2) +
  #   # stat_summary(fun.y=median, geom="line", size = 2, color = "blue") +
  #   labs(title = paste(stimulus, '_', test_gene1_name)) + guides(color = FALSE) +
  #   scale_color_discrete()+theme_bw()
  #    
  
  # print(plot_points)
  print(plot_lines)
  # print(plot_few)
  # print(plot_density)
  # print(plot_lines_sd)
  
  # toRet[['plot1']] = plot_points
  toRet[['plot2']] = plot_lines
  # toRet[['plot3']] = plot_few
  # toRet[['plot4']] = plot_density
  # toRet[['plot5']] = plot_lines_sd
  
  if(scale != 'none'){
    toRet[['scaled_gene_exp']] = test_gene1$exp
  }
  else{
    toRet[['scaled_gene_exp']] = 'Scaling was set to FALSE'
  }
  return(toRet)
}

getMeanMedVar_NoPlots = function(macro, metadata, reconstructed_data, gene, timepoints, scale, omitzero){
  require('matrixStats')
  toRet = list()
  
  gene_ind = grep(pattern = gene, x = colnames(reconstructed_data[['reconstructed_trajectories']]))
  #print(gene_ind)
  #print(colnames(reconstructed_data[['reconstructed_trajectories']])[gene_ind])
  if(length(gene_ind)==0){
    print('Gene not found in simulated data')
    return(NULL)
  }
  else{
    gene_ind = gene_ind[1] # uses first instance of desired gene
  }
  test_gene1 = reconstructed_data[['reconstructed_trajectories']][,c(gene_ind,ncol(reconstructed_data$reconstructed_trajectories)-1, ncol(reconstructed_data$reconstructed_trajectories))]
  #View(test_gene1)
  colnames(test_gene1)[1] = 'exp'
  if(scale){
    test_gene1$exp = (test_gene1$exp - min(test_gene1$exp))/(max(test_gene1$exp)-min(test_gene1$exp))
  }
  
  truedata = as.data.frame(macro@assays$ISnorm@data)
  gene_ind_true = grep(pattern = gene, x = rownames(truedata))
  if(length(gene_ind_true)==0){
    print('Gene not found in true data')
    return(NULL)
  }
  else{
    gene_ind_true = gene_ind_true[1] # uses first instance of desired gene
  }
  # print('these are cols')
  # print(colnames(truedata))
  # 
  # print('these are rows')
  # print(rownames(truedata))
  truedata = truedata[gene_ind_true,]
  
  if(omitzero){
    #truedata = truedata[,truedata != 0]
    truedata = truedata[,-(which(colSums(truedata)==0))]
  }
  
  if(is.null(dim(truedata)) || 0 %in% dim(truedata)){
    print(paste("the following gene has been omitted due to insufficient data after filtering:", gene))
    return(NULL)
  }
  
  if(scale){
    truedata[1,] = (truedata[1,]-min(truedata[1,]))/(max(truedata[1,])- min(truedata[1,]))
  }
  
  #View(truedata)
  #View(test_gene1)
  sim_means = c()
  sim_meds = c()
  sim_var = c()
  data_means = c()
  data_var = c()
  data_meds = c()
  #boxplots = list()
  
  for(time in timepoints){
    timeptstr = paste(time, 'hr', sep = "")
    #print(paste('Checking time', timeptstr))
    
    min_time = test_gene1$time[which.min(abs(time-0.5-test_gene1$time))]
    max_time = test_gene1$time[which.min(abs(time+0.5-test_gene1$time))]
    times = unique(test_gene1$time[test_gene1$time > min_time & test_gene1$time < max_time])
    #print(times)
    
    relevant_cells = rownames(metadata)[metadata$timept_num == time] # Get all cells that were measured at given time point 
    #print(length(relevant_cells))
    relevant_cells =  colnames(truedata) %in% relevant_cells
    
    truecells = cbind(t(truedata[1,relevant_cells]), 'True_Data')
    simcells = cbind(test_gene1$exp[test_gene1$time %in% times], 'Simulated_Data')
    plottabledf = as.data.frame(rbind(truecells, simcells))
    colnames(plottabledf) = c('value', 'group')
    plottabledf$value = as.numeric(plottabledf$value)
    plottabledf$group = as.factor(plottabledf$group)
    
    #View(plottabledf)
    #print(str(plottabledf))
    #boxplots[[timeptstr]] = ggplot(plottabledf, aes(x = group, y = value, fill = group)) + geom_boxplot() + ggtitle(paste('Comparison of simulated \n and true distributions of', gene, 'at', timeptstr)) + ylab(ifelse(scale, "scaled value", "unscaled value")) + geom_jitter(alpha=0.2)
    
    #View(as.matrix(truedata[1, relevant_cells %in% colnames(truedata)]))
    #datamean_cur = as.numeric(rowMeans(as.matrix(truedata[1, relevant_cells]))) # Get mean of true data for all relevant cells
    datamean_cur = as.numeric(mean(plottabledf$value[plottabledf$group == 'True_Data']))
    datamed_cur = as.numeric(median(plottabledf$value[plottabledf$group == 'True_Data']))
    datavar_cur = as.numeric(var(plottabledf$value[plottabledf$group == 'True_Data']))
    #datavar_cur = as.numeric(rowVars(as.matrix(truedata[1, relevant_cells]))) # Get variance of true data for all relevant cells
    data_means = c(data_means, datamean_cur)
    data_meds = c(data_meds, datamed_cur)
    data_var = c(data_var, datavar_cur)
    #print(data_means)
    
    #sim_means = c(sim_means, mean(test_gene1$exp[test_gene1$time == time]))
    sim_means = c(sim_means, as.numeric(mean(plottabledf$value[plottabledf$group == 'Simulated_Data'])))
    sim_meds = c(sim_meds, as.numeric(median(plottabledf$value[plottabledf$group == 'Simulated_Data'])))
    sim_var = c(sim_var, as.numeric(var(plottabledf$value[plottabledf$group == 'Simulated_Data'])))
    #print(sim_means)
    #sim_var = c(sim_var, var(test_gene1$exp[test_gene1$time == time]))
  }
  
  ret_meds = c(data_meds, sim_meds)
  ret_vars = c(data_var, sim_var)
  ret_means = c(data_means, sim_means)
  #ret = cbind(paste(timepoints, 'hrs', sep = ''), data_means)
  #ret = cbind(ret, data_var)
  #ret = cbind(ret, sim_means)
  #ret = cbind(ret, sim_var)
  #colnames(ret) = c('Timepoint', paste('True mean', gene, sep = '-'), paste('True variance', gene, sep = '-'), paste('Simulated mean', gene, sep = '-'), paste('Simulated variances', gene, sep = '-'))
  
  #toRet[['stats']] = ret
  #toRet[['boxplots']] = boxplots
  #toRet[['gene']] = gene
  
  toRet[['means']] = ret_means
  toRet[['medians']] = ret_meds
  toRet[['variances']] = ret_vars
  
  
  return(toRet)
}
