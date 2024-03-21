#' UG_methods
#'
#' UG_methods is the function that use indirect graph algorithm to get the connection between nodes
#'
#' @param mem matrix of data, gene as rows and cell types as cols
#' @param param numeric, the parameter of UG methods, ranging from 0 to 1
#' @param plot logic
#' @param weight logic
#' @param ugMethod character, determine the method used, usually cmi2ni. Other alternatives including ns/GeneNet/glasso/pcacmi/bayesianglasso
#'
#' @return matrix, return the UG method calculation results
#' @export
#'
#' @examples
#' test=ma_test[names(ma_test_glist)[1:3000],]
#' UG_methods(test,param=0.2)
UG_methods <- function(mem, param = 0.2, plot = FALSE, weight = FALSE, ugMethod = 'cmi2ni'){
  if (ugMethod == 'GeneNet') {
    ug <- network.GeneNet(expr.data = (mem), fdr = param)
  } else if (ugMethod == 'ns') {
    ug <- network.ns(expr.data = (mem), alpha = param)
  } else if (ugMethod == 'glasso') {
    ug <- network.glasso(expr.data = (mem), lambda = param)
  } else if (ugMethod == 'glassosf') {
    ug <- network.glassosf(expr.data = (mem), alpha = param)
  } else if (ugMethod == 'pcacmi') {
    ug <- network.pcacmi(expr.data = (mem), lambda = param)
  } else if (ugMethod == 'cmi2ni') {
    ug <- network.cmi2ni(expr.data = (mem), lambda = param)
  } else if (ugMethod == 'space') {
    ug <- network.space(expr.data = (mem), alpha = param)
  } else if (ugMethod == 'bayesianglasso') {
    ug <- network.bayesianglasso(expr.data = (mem), prob = param)
  }
  
  if(plot) ugPlot(ug,main = paste0(ugMethod,', para = ', param),
                  weight = weight)
  return(ug)
}


#' DG_grid
#'
#' DG_grid is the function that considers all the parameters in "params" to get a group of net structures based on ugMethod and dagMethod
#'
#' @param mem matrix of data, gene as rows and cell types as cols
#' @param params numeric, the parameters of UG methods, ranging from 0 to 1
#' @param whiteList data.frame, the list of arcs that must exist in the network
#' @param blackList data.frame, the list of arcs that can't exist in the network
#' @param root character, the root of the graph, which can't have any node pointing to.
#' @param ugMethod character, determine the UG_method used, usually cmi2ni. Other alternatives including ns/GeneNet/glasso/pcacmi/bayesianglasso
#' @param dagMethod character, determine the method that calculate the arc direction, usually hc
#' @param ncores integer,the cores that used to speed up the calculation
#'
#' @return A list of net structures of different parameters
#' @export
#'
#' @examples
#' library(foreach)
#' library(bnlearn)
#' test <- ma_test[names(ma_test_glist)[1:3000],]
#' net_struc_ori <- DG_grid(test,params=seq(0.08,0.23,0.01),root='HSC')
#' graphviz.plot(net_struc_ori[[1]])
DG_grid <- function(mem, params = c(1:100)/100, whiteList = NULL, blackList = NULL, root = NULL,
                    ugMethod = 'cmi2ni', dagMethod = 'hc', ncores = 1){
  doParallel::registerDoParallel(cores = ncores)
  dags <- foreach::foreach(param = params, .combine = c) %dopar% {
    tmp <- list()
    tmp[[1]] <- DG_methods(mem = (mem), param = param, root = root,
                           whiteList = whiteList, blackList = blackList,
                           ugMethod = ugMethod, dagMethod = dagMethod)
    tmp
  }
  names(dags) <- paste0('param_',params)
  dags
}

#' DG_methods
#'
#' DG_methods is the function that use one parameter to get a directed net struture based on ugMethod and dagMethod
#'
#' @param mem matrix of data, gene as rows and cell types as cols
#' @param param numeric, the parameters of UG methods, ranging from 0 to 1
#' @param root character, the root of the graph, which can't have any node pointing to.
#' @param whiteList data.frame, the list of arcs that must exist in the network
#' @param blackList data.frame, the list of arcs that can't exist in the network
#' @param plot logic, if Ture, the net structure of this parameter will be ploted
#' @param ugMethod character, determine the UG_method used, usually cmi2ni. Other alternatives including ns/GeneNet/glasso/pcacmi/bayesianglasso
#' @param dagMethod character, determine the method that calculate the arc direction, usually hc
#'
#' @return The net structure learnt under this parameter and method
#' @export
#'
#' @examples
#' library(foreach)
#' library(bnlearn)
#' test <- ma_test[names(ma_test_glist)[1:3000],]
#' result <- DG_methods(test,param=0.08,root='HSC')
#' graphviz.plot(result)
DG_methods <- function(mem, param = 0.2, root = NULL, whiteList = NULL, blackList = NULL, plot = FALSE,
                       ugMethod = 'cmi2ni', dagMethod = 'hc'){
  "Edges in whiteList and blackList should be in the form of node1~node2"
  ug <- UG_methods((mem), param = param, ugMethod = ugMethod)
  nodes <- colnames((mem))
  blackList <- genBList(ug, nodes = nodes, root = root, blackList = blackList)
  dag <- get(dagMethod)(as.data.frame(mem), blacklist = blackList, whitelist = whiteList)
  main = paste0('method = ',ugMethod,', ng = ', nrow((mem)),', param = ',param)
  if (plot) bnlearn::graphviz.plot(dag, main = main,shape = "circle")
  dag
}


#' DG_smpl
#'
#' DG_smpl is the function that considers all the parameters in "params" and cell samplings to get a group of net structures based on ugMethod and dagMethod. Compared with DG_grid, DG_smpl is suitable for single cell data in SeuratObject format
#'
#' @param dat Seurat object for input
#' @param frac numeric, the fraction of cells for every sampling
#' @param N_smpl numeric, the total times of sampling
#' @param params numeric, the parameters used for UG_method
#' @param whiteList data.frame, the list of arcs that must exist in the network
#' @param blackList data.frame, the list of arcs that can't exist in the network
#' @param root character, the root of the graph, which can't have any node pointing to.
#' @param ugMethod character, determine the UG_method used, usually cmi2ni. Other alternatives including ns/GeneNet/glasso/pcacmi/bayesianglasso
#' @param dagMethod character, determine the method that calculate the arc direction, usually hc
#' @param ncores integer,the cores that used to speed up the calculation
#'
#' @return A list of net structures of different parameters and different sampling
#' @export
#'
#' @examples
#' library(foreach)
#' library(bnlearn)
#' library(Seurat)
#' dat <- CreateSeuratObject(dat_fibro)
#' dat$celltype <- dat_fibro_meta
#' Idents(dat)=dat$celltype
#' a <- DG_smpl(dat[fibro_glist,], frac = 0.1, N_smpl = 10, params = seq(0.001,0.01,0.001),whiteList = NULL, blackList = NULL, root = 'MEF',ugMethod = 'cmi2ni', dagMethod = 'hc', ncores = 64)
#' graphviz.plot(a[[1]][[1]])
DG_smpl <- function (dat, frac = 0.3, N_smpl = 10, params = c(1:100)/100, 
                     whiteList = NULL, blackList = NULL, root = NULL, ugMethod = "cmi2ni", 
                     dagMethod = "hc", ncores = 64,seed=Sys.time()) 
{
  set.seed(as.numeric(seed))
  meta <- dat@meta.data %>% tibble::rownames_to_column("CellID")
  ctypes = names(table(meta$celltype))
  random_number = runif(N_smpl, 1, as.numeric(seed))
  
  dags.smpl <- foreach::foreach(i_smpl = 1:N_smpl, .combine = c) %do% 
    {
      set.seed(random_number[i_smpl])
      meta.smpl <- meta %>% dplyr::group_by(celltype) %>% 
        dplyr::sample_frac(frac)
      dat.smpl <- subset(dat, cells = meta.smpl$CellID)
      meta_tmp = dat.smpl$celltype
      tmp = data.frame(dat.smpl@assays$RNA@data)
      mem.smpl <- data.frame(matrix(0, nrow = nrow(tmp), 
                                    ncol = length(ctypes)))
      colnames(mem.smpl) <- ctypes
      rownames(mem.smpl) <- rownames(tmp)
      for (c_t in ctypes) {
        mem.smpl[c_t] <- tmp[, meta_tmp == c_t] %>% 
          rowMeans
      }
      dags <- list()
      # rm(c)
      dags[[1]] <- DG_grid(mem.smpl, params = params, 
                           whiteList = whiteList, blackList = blackList, 
                           root = root, ugMethod = ugMethod, dagMethod = dagMethod, 
                           ncores = ncores)
      # dags[[2]] <- meta.smpl
      # dags[[3]] <- seed
      dags
    }
  dags.smpl
}

#' BNLearning
#'
#' BNLearning is the function that considers all the parameters in "params" and cell samplings to get a group of net structures based on ugMethod and dagMethod. Compared with DG_grid, DG_smpl is suitable for single cell data in SeuratObject format
#'
#' @param dat Seurat object for input
#' @param frac numeric, the fraction of cells for every sampling
#' @param N_smpl numeric, the total times of sampling
#' @param params numeric, the parameters used for UG_method
#' @param whiteList data.frame, the list of arcs that must exist in the network
#' @param blackList data.frame, the list of arcs that can't exist in the network
#' @param root character, the root of the graph, which can't have any node pointing to.
#' @param ugMethod character, determine the UG_method used, usually cmi2ni. Other alternatives including ns/GeneNet/glasso/pcacmi/bayesianglasso
#' @param dagMethod character, determine the method that calculate the arc direction, usually hc
#' @param ncores integer,the cores that used to speed up the calculation
#'
#' @return A list of net structures of different parameters and different sampling
#' @export
#'
#' @examples
#' library(foreach)
#' library(bnlearn)
#' library(Seurat)
#' dat <- CreateSeuratObject(dat_fibro)
#' dat$celltype <- dat_fibro_meta
#' Idents(dat)=dat$celltype
#' a <- BNLearning(dat[fibro_glist,], frac = 0.1, N_smpl = 10, params = seq(0.001,0.01,0.001),whiteList = NULL, blackList = NULL, root = 'MEF',ugMethod = 'cmi2ni', dagMethod = 'hc', ncores = 1,mode='single_cell')
#' graphviz.plot(a[[1]][[1]])
BNLearning <- function(dat, frac = 0.3, N_smpl = 10, params = c(1:100)/100,
                       whiteList = NULL, blackList = NULL, root = NULL,
                       ugMethod = 'cmi2ni', dagMethod = 'hc', ncores = 1,
                       mode = 'single_cell'){
  if (mode == 'single_cell') {
    BNLearn_result <- DG_smpl(dat, frac , N_smpl , params ,
                              whiteList , blackList , root ,
                              ugMethod , dagMethod , ncores)
  } else if (mode == 'bulk') {
    BNLearn_result <- DG_grid(dat, params , whiteList , blackList , root ,
                              ugMethod , dagMethod , ncores )
  } else {BNLearn_result=NULL; print('Wrong mode!')}
  BNLearn_result
}


#' Title bootstrap_index
#'
#' bootstrap_index: create sampling index for run_diffBN
#'
#' @param meta character, the cell type information of single cell data
#' @param bootstrap_times numeric, times of sampling, should be some with "bootstrap_time" in "run_diffBN"
#' @param ratio double, the fraction of cell that are involved in every sample
#'
#' @return list of index for diffBN
#' @export
#'
#' @examples
#' index <- bootstrap_index(dat_fibro_meta,bootstrap_times=20,ratio=0.4)
bootstrap_index <- function(meta,bootstrap_times,ratio){
  celltype=names(table(meta))
  index=list(NULL)
  for (j in 1:bootstrap_times)
  {
    index[[j]]=NA
    for (i in 1:length(table(meta)))
    {
      tmp=which(meta==celltype[i])
      index[[j]]=c(index[[j]],sample(tmp,length(tmp)*ratio,replace = F))
    }
    index[[j]]=index[[j]][-1]
  }
  
  return(index)
}

#' Title combineDAGsmpl
#'
#' combineDAGsmpl is a function that combine all the net structures in the dags.smpl according to a filter of Emin and Emax to get the comprehensize network.
#'
#' @param dags.smpl list of net structures, the result of DG_smpl
#' @param Emin integer, the minimal arc number. Any net that have arc number lower that Emin will be filtered
#' @param Emax integer, the maximal arc number. Any net that have arc number higher that Emax will be filtered
#' @param ncores integer, number of cores used to speed up the calculation
#'
#' @return information of the summarized net structure
#' @export
#'
#' @examples
#' library(foreach)
#' library(bnlearn)
#' library(Seurat)
#' dat <- CreateSeuratObject(dat_fibro)
#' dat$celltype <- dat_fibro_meta
#' Idents(dat)=dat$celltype
#' a <- DG_smpl(dat[fibro_glist,], frac = 0.1, N_smpl = 10, params = seq(0.001,0.01,0.001),whiteList = NULL, blackList = NULL, root = 'MEF',ugMethod = 'cmi2ni', dagMethod = 'hc', ncores = 64)
#' b <- combineDAGsmpl(a,Emin=1,Emax=8)
combineDAGsmpl <- function(dags.smpl, Emin = NULL, Emax = NULL, ncores = 64){
  Efreqs.tmp <- foreach::foreach(i_smpl = 1:length(dags.smpl), .combine = rbind) %dopar% {
    combineDAGs(dags.smpl[[i_smpl]], Emin = Emin, Emax = Emax)
  }
  Efreqs.tmp %>% dplyr::group_by(from, to) %>% dplyr::summarise(freq = sum(freq)) %>% as.data.frame %>%
    dplyr::mutate(edge = paste0(from,'~',to)) %>% tibble::remove_rownames() %>% tibble::column_to_rownames('edge')
}




#' Title combineDAGs
#'
#' combineDAGs is a function that combine all the net structures in the dags according to a filter of Emin and Emax to get the comprehensize network.
#'
#' @param dags list of net structures, the result of DG_grid
#' @param Emin integer, the minimal arc number. Any net that have arc number lower that Emin will be filtered
#' @param Emax integer, the maximal arc number. Any net that have arc number higher that Emax will be filtered
#'
#' @return information of the summarized net structure
#' @export
#'
#' @examples
#' library(foreach)
#' library(bnlearn)
#' test <- ma_test[names(ma_test_glist)[1:3000],]
#' net_struc_ori <- DG_grid(test,params=seq(0.08,0.23,0.01),root='HSC')
#' net_struc <- combineDAGs(net_struc_ori,Emin=20,Emax=30)
combineDAGs <- function(dags, Emin = NULL, Emax = NULL){
  Efreqs <- foreach::foreach(i = 1:length(dags), .combine = rbind) %do% {
    bool1 <- ifelse(is.null(Emin), TRUE, bnlearn::narcs(dags[[i]]) >= Emin)
    bool2 <- ifelse(is.null(Emax), TRUE, bnlearn::narcs(dags[[i]]) <= Emax)
    bool <- bool1 && bool2
    # bool <- (bnlearn::narcs(dags[[i]]) >= Emin) && (bnlearn::narcs(dags[[i]]) <= Emax)
    if (isTRUE(bool)) {
      Efreq <- bnlearn::arcs(dags[[i]]) %>% as.data.frame %>% dplyr::mutate(freq = 1)
      Efreq
    } else NULL
  }
  Efreqs %>% dplyr::group_by(from, to) %>% dplyr::summarise(freq = sum(freq)) %>% as.data.frame %>%
    dplyr::mutate(edge = paste0(from,'~',to)) %>% tibble::remove_rownames() %>% tibble::column_to_rownames('edge')
}

#' Title trimDAG
#'
#' trimDAG is a function that modify netword structure according to bn.fit results
#'
#' @param dat_tmp matrix of gene expression based on different cell types
#' @param e A bn.fit net structure
#' @param min_arc integer, the minimal number of arcs pointing to a node
#' @param max_arc integer, the maximal number of arcs pointing to a node
#' @param threshold_value numeric, determine the filter threshold with which arcs are filtered according to bn.fit results
#'
#' @return A new bn.fit net structure after filter
#' @export
#'
#' @examples
#' library(bnlearn)
#' e <-  empty.graph(colnames(ma_test))
#' net_struc <- rmCyc(net_struc) # net_struc is the result of combineDAGs
#' strc1 <- (net_struc %>% rmCyc)[,1:2]
#' strc1[,1]=as.character(strc1[,1])
#' strc1[,2]=as.character(strc1[,2])
#' arcs(e) <-  strc1
#' trimDAG(ma_test,e,2,4,1)
trimDAG <- function(dat_tmp,e,min_arc = 2,max_arc = 4,threshold_value = 0.9){
  
  ref = bnlearn::bn.fit(e,dat_tmp)
  celltype = names(ref)
  arc_modified = bnlearn::arcs(ref)
  
  # cutoff & min
  for (i in celltype)
  {
    cut_off = sum(abs(unlist(ref[[i]][4])[-1]))*(1-threshold_value)
    if (length(unlist(ref[[i]][4])[-1]) == 0) {i = i}
    else if (length(unlist(ref[[i]][4]))<(min_arc+2)) {i = i}
    else {
      for (j in 1:length(unlist(ref[[i]][4])[-1]))
      {if (abs(unlist(ref[[i]][4])[-1])[j]<cut_off)
      {
        tmp = strsplit(names(unlist(ref[[i]][4])[-1])[j],split = 'coefficients.')[[1]][2]
        target = intersect(which(arc_modified[,'to'] == i),which(arc_modified[,'from'] == tmp))
        arc_modified = arc_modified[-target,]
      }
      }
    }
  }
  bnlearn::arcs(e) = arc_modified
  ref = bnlearn::bn.fit(e,dat_tmp)
  
  # max
  for (i in celltype)
  {
    if (length(unlist(ref[[i]][4])[-1])>max_arc)
    {
      cut_off = as.numeric(sort(abs(unlist(ref[[i]][4])[-1]),decreasing = T)[max_arc])
      for (j in 1:length(unlist(ref[[i]][4])[-1]))
      {if (abs(unlist(ref[[i]][4])[-1])[j]<cut_off)
      {
        tmp = strsplit(names(unlist(ref[[i]][4])[-1])[j],split = 'coefficients.')[[1]][2]
        target = intersect(which(arc_modified[,'to'] == i),which(arc_modified[,'from'] == tmp))
        arc_modified = arc_modified[-target,]
      }
      }
    }
  }
  bnlearn::arcs(e) = arc_modified
  bnlearn::graphviz.plot(e, shape = "ellipse")
  return(e)
}


#' Title arc_modify
#'
#' arc_modify is a function that modify netword structure according to bn.fit results
#'
#' @param dat_tmp matrix of gene expression based on different cell types
#' @param e A bn.fit net structure
#' @param min_arc integer, the minimal number of arcs pointing to a node
#' @param max_arc integer, the maximal number of arcs pointing to a node
#' @param threshold_value numeric, determine the filter threshold with which arcs are filtered according to bn.fit results
#'
#' @return A new bn.fit net structure after filter
#' @export
#'
#' @examples
#' library(bnlearn)
#' e <-  empty.graph(colnames(ma_test))
#' strc1 <- (net_struc %>% rmCyc)[,1:2]
#' strc1[,1]=as.character(strc1[,1])
#' strc1[,2]=as.character(strc1[,2])
#' arcs(e) <-  strc1
#' arc_modify(ma_test,e,2,3,1)
arc_modify <- function(dat_tmp,e,min_arc = 2,max_arc = 4,threshold_value = 0.9){
  
  ref = bnlearn::bn.fit(e,dat_tmp)
  celltype = names(ref)
  arc_modified = bnlearn::arcs(ref)
  
  # cutoff & min
  for (i in celltype)
  {
    cut_off = sum(abs(unlist(ref[[i]][4])[-1]))*(1-threshold_value)
    if (length(unlist(ref[[i]][4])[-1]) == 0) {i = i}
    else if (length(unlist(ref[[i]][4]))<(min_arc+2)) {i = i}
    else {
      for (j in 1:length(unlist(ref[[i]][4])[-1]))
      {if (abs(unlist(ref[[i]][4])[-1])[j]<cut_off)
      {
        tmp = strsplit(names(unlist(ref[[i]][4])[-1])[j],split = 'coefficients.')[[1]][2]
        target = intersect(which(arc_modified[,'to'] == i),which(arc_modified[,'from'] == tmp))
        arc_modified = arc_modified[-target,]
      }
      }
    }
  }
  bnlearn::arcs(e) = arc_modified
  ref = bnlearn::bn.fit(e,dat_tmp)
  
  # max
  for (i in celltype)
  {
    if (length(unlist(ref[[i]][4])[-1])>max_arc)
    {
      cut_off = as.numeric(sort(abs(unlist(ref[[i]][4])[-1]),decreasing = T)[max_arc])
      for (j in 1:length(unlist(ref[[i]][4])[-1]))
      {if (abs(unlist(ref[[i]][4])[-1])[j]<cut_off)
      {
        tmp = strsplit(names(unlist(ref[[i]][4])[-1])[j],split = 'coefficients.')[[1]][2]
        target = intersect(which(arc_modified[,'to'] == i),which(arc_modified[,'from'] == tmp))
        arc_modified = arc_modified[-target,]
      }
      }
    }
  }
  bnlearn::arcs(e) = arc_modified
  bnlearn::graphviz.plot(e, shape = "ellipse")
  return(e)
}

#' Title rmCycL
#'
#' rmCycL: remove the circles in the network according to the hierarchical structure
#'
#' @param dS Summarized net structure
#'
#' @return Summarized net structure, but the circles are removed
#' @export
#'
#' @examples
#' rmCycL(net_struc) # net_struc is the result of combineDAGs
rmCycL <- function(dS){
  Es <- dS[c('from','to')]
  Es[,1]=as.character(Es[,1])
  Es[,2]=as.character(Es[,2])
  # gp <- Es %>% as.matrix %>% graph::ftM2graphNEL
  # lgp <- Rgraphviz::layoutGraph(gp)
  # Ys <- graph::nodeRenderInfo(lgp)$nodeY
  Ys <- Es %>% as.matrix %>% (graph::ftM2graphNEL) %>%
    (Rgraphviz::layoutGraph) %>% (graph::nodeRenderInfo) %>% .$nodeY %>% sort %>% rev
  coords <- Es %>% dplyr::mutate(from = Ys[from], to = Ys[to])
  dS.fltr <- dS[coords$from > coords$to,]
  dS.fltr
}

#' Title
#'
#' @param dS list
#'
#' @return list
#' @export
#'
#' @examples
#' print('ok')
rmCyc2 <- function(dS){
  dS <- dS %>% df2mat
  dS[dS-t(dS)<0] = 0
  dS %>% mat2df
}

#' Title rmCyc
#'
#' rmCyc: remove the circles in the network according to the hierarchical structure and net matrix
#'
#' @param dS Summarized net structure
#'
#' @return Summarized net structure, but the circles are removed
#' @export
#'
#' @examples
#' rmCyc(net_struc) # net_struc is the result of combineDAGs
rmCyc <- function(dS){
  tmp <- dS %>% rmCyc2
  dS.fltr <- tmp %>% rmCycL
  edge_rev <- dS.fltr %>% dplyr::mutate(edge_rev = paste0(to,'~',from)) %>% .$edge_rev
  edge_tune <- setdiff(rownames(dS), rownames(tmp)) %>% setdiff(., edge_rev) %>% union(., rownames(dS.fltr))
  dS.tune <- dS[edge_tune,]
  dS.tune
}


#' Title
#'
#' @param N matrix
#'
#' @return matrix
#' @export
#'
#' @examples
#' print('ok')
norm_mat <- function(N){
  M <- N
  for(c in 1:nrow(M)){
    for(d in 1:ncol(M)){
      if(c<d){
        if(M[c,d] + M[d,c] > 0 ){
          if(M[c,d] > M[d,c]){
            M[c,d] <- 1
            M[d,c] <- 0
          }else{
            M[c,d] <- 0
            M[d,c] <- 1
          }
        }
      }
    }
  }
  return(M)
}

#' Title
#'
#' @param fromSet character
#' @param toSet character
#' @param sep character
#'
#' @return character
#' @export
#'
#' @examples
#' print('ok')
setEdges <- function(fromSet, toSet, sep = '~'){
  edges <- paste(rep(fromSet, each = length(toSet)), toSet, sep = sep)
  selfs <- paste(fromSet, fromSet, sep = sep)
  return(setdiff(edges,selfs))
}


#' gem2mem
#'
#' gem2mem is a function that can summary a single cell expression matrix to a cell type expression matrix using the metadata and method given
#' Usually, FUN='mean'
#'
#' @param gem data.frame of single cell expression
#' @param meta character, cell type information of single cell expression data
#' @param FUN character, can be 'mean' or 'median'
#'
#' @return data.frame
#' @export
#'
#' @examples
#' fibro_celltype <- gem2mem(dat_fibro,dat_fibro_meta,FUN='mean')
#' View(fibro_celltype)
gem2mem <- function(gem=NULL, meta=NULL, FUN=c('mean','median')){
  if (is.null(meta)) {
    warning('Celltypes information is missing, MEM is identical to input GEM ...')
    mem <- gem
  } else{
    celltype <- names(table(meta))
    mem = data.frame(matrix(0,nrow=nrow(gem),ncol=length(celltype)))
    colnames(mem) = celltype
    rownames(mem) = rownames(gem)
    for (c in celltype)
    {
      if (FUN=='mean') {mem[,c]=rowMeans(gem[,meta==c])}
      else if (FUN=='median') {mem[,c]=matrixStats::rowMedians(gem[,meta==c])}
    }
  }
  return(mem)
}

#' Title
#'
#' @param nEdges interger
#' @param Emin interger
#' @param Emax interger
#'
#' @return interger
#' @export
#'
#' @examples
#' print('ok')
filterParams <- function(nEdges, Emin = min(nEdges), Emax = max(nEdges)){
  mask <- (nEdges >= Emin) & (nEdges <= Emax)
  params <- names(nEdges)[mask] %>% as.numeric
  params
}

#' Title
#'
#' @param ug data.frame
#' @param nodes character
#' @param root character
#' @param blackList data.frame
#'
#' @return data.frame
#' @export
#'
#' @examples
#' print('he')
genBList <- function(ug, nodes, root = NULL, blackList = NULL){
  edgeAll <- setEdges(nodes,nodes,sep = "~")
  edgeCor <- paste0(c(ug$node1,ug$node2),'~',c(ug$node2,ug$node1))
  blacklist <- setdiff(edgeAll, edgeCor) %>%
    union(., setEdges(nodes, root)) %>%
    union(., blackList) %>%
    as.data.frame %>% rlang::set_names('edge') %>%
    tidyr::separate(edge, c('from','to'), sep = '~')
  blacklist
}

#' Title
#'
#' @param nEdges interger
#' @param Emin interger
#' @param Emax interger
#'
#' @return interger
#' @export
#'
#' @examples
#' print('ok')
filterParams <- function(nEdges, Emin = min(nEdges), Emax = max(nEdges)){
  mask <- (nEdges >= Emin) & (nEdges <= Emax)
  params <- names(nEdges)[mask] %>% as.numeric
  params
}

#' Title
#'
#' @param ... character
#'
#' @return character
#' @export
#'
#' @examples
#' print('ok')
alter <- function(...){
  vec <- unique(c(rbind(...)))
  return(vec)
}

#' Title
#'
#' @param dags list
#'
#' @return character
#' @export
#'
#' @examples
#' print('ok')
getNEdges <- function(dags){
  nEdges <- sapply(1:length(dags),
                   function(x) bnlearn::narcs(dags[[x]])) %>% rlang::set_names(names(dags))
  nEdges
}

#' Title df2mat
#'
#' df2mat: transform a summarized net structure to matrix
#'
#' @param DAG Summarized net structure
#' @param ctypes character, the cell typs used
#'
#' @return A matrix that is corresponding to the summarized net structure
#' @export
#'
#' @examples
#' library(foreach)
#' library(bnlearn)
#' test <- ma_test[names(ma_test_glist)[1:3000],]
#' net_struc_ori <- DG_grid(test,params=seq(0.08,0.23,0.01),root='HSC')
#' net_struc <- combineDAGs(net_struc_ori,Emin=20,Emax=30)
#' df2mat(net_struc)
df2mat <- function(DAG, ctypes = NULL){
  DAG.dm <- tidyr::spread(DAG, key = to, value = freq) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("from") %>%
    dplyr::mutate(dplyr::across(tidyselect::everything(), .fns = ~replace_na(.,0)))
  
  if (is.null(ctypes)) ctypes <- union(rownames(DAG.dm), colnames(DAG.dm))
  DAG.dm[setdiff(ctypes,rownames(DAG.dm)),] <- 0
  DAG.dm[,setdiff(ctypes,colnames(DAG.dm))] <- 0
  DAG.dm <- DAG.dm[ctypes,ctypes]
  DAG.dm
}

#' Title mat2df
#'
#' mat2df: transform a matrix to a summarized net structure
#'
#' @param DAG.dm A matrix that represent a net structure
#'
#' @return Summarized net structure
#' @export
#'
#' @examples
#' library(foreach)
#' library(bnlearn)
#' test <- ma_test[names(ma_test_glist)[1:3000],]
#' net_struc_ori <- DG_grid(test,params=seq(0.08,0.23,0.01),root='HSC')
#' net_struc <- combineDAGs(net_struc_ori,Emin=20,Emax=30)
#' net_struc_mat <- df2mat(net_struc)
#' net_struc <- mat2df(net_struc_mat)
mat2df <- function(DAG.dm){
  DAG <- DAG.dm  %>% t %>% as.table %>% as.data.frame %>%
    rlang::set_names(c("to","from","freq")) %>% subset(freq != 0) %>%
    .[,c("from","to","freq")] %>%
    dplyr::mutate(from = as.character(from), to = as.character(to), edge = paste0(from,'~',to)) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames('edge')
  DAG
}
