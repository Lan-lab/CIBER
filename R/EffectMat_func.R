#' Title calculate_diffBN
#'
#' When function "run_diffBN" is done, the ouput of "run_diffBN" can be the input of this function to produce diffBN results of other modes
#'
#' @param diffBN list of diffBN raw result
#' @param diffBN_mode character, determine the diffBN function, usually "mean"
#'
#' @return diffBN results
#' @export
#'
#' @examples
#' calculate_diffBN(test_result,diffBN_mode='mean') # test_result is the result of function "run_diffBN"
calculate_diffBN <- function(diffBN,diffBN_mode='mean')
{
  diffBN_result <- list(NULL)
  bootstrap_time = diffBN$bootstrap_time
  permutation_time = diffBN$permutation_time
  if (diffBN_mode == 'mean'){
    k=1
    for (i in 1:(bootstrap_time))
    {
      for (j in 1:(permutation_time))
      {
        diffBN_result[[k]] = get_diffBN(scm.inter = diffBN$raw[[k]], scm.ref = diffBN$ref[[i]], mode = diffBN_mode)
        k=k+1
      }
    }
    names(diffBN_result) = names(diffBN$raw)
  } else if (diffBN_mode == 'OT'){
    n_edge <- nrow(diffBN$ref[[1]])
    n_features <- ncol(diffBN$raw[[1]])
    diffBN_ref <- array(unlist(diffBN$ref),dim = c(bootstrap_time,n_edge))
    diffBN_raw <- array(unlist(diffBN$raw),dim = c(bootstrap_time*permutation_time,n_features,n_edge))
    # row=gene col=edge
    diffBN_result <- foreach(i = 1:n_features, .combine = rbind) %dopar%{
      foreach(j = 1:n_edge, .combine = c) %dopar%{
        wasserstein1d(diffBN_ref[,j],diffBN_raw[,i,j])
      }
    }
    colnames(diffBN_result) <- rownames(diffBN$ref[[1]])
    rownames(diffBN_result) <- colnames(diffBN$raw[[1]])
  } else{
    stop("Invalid diffBN calculation mode!")
  }
  return(as.data.frame(diffBN_result))
}

#' EffectMatrix
#'
#' EffectMatrix is a function that calculate diffBN results based on a given expression data
#'
#' @param net_struc bn.fit structure
#' @param data data.frame of the expression data
#' @param meta character of metadata corresponding to expression data
#' @param index list of sampling results
#' @param permutation_time integer, times of gene permutation, usually 20-50
#' @param bootstrap_time integer, times of cell sampling, usually 20-50
#' @param diffBN_mode character, determine the diffBN function, usually "mean"
#'
#' @return diffBN results
#' @export
#'
#' @examples
#' library(foreach)
#' # not run: test_result <- EffectMatrix(net_struc,dat_fibro,dat_fibro_meta,index,permutation_time=20,bootstrap_time=20,diffBN_mode='mean',mode='single_cell)
EffectMatrix <- function(net_struc,data,
                         meta=NULL,index=NULL,permutation_time=20,bootstrap_time=20,
                         diffBN_mode='mean',mode='single_cell',ncores=1)
{
  doParallel::registerDoParallel(ncores)
  if (mode == 'single_cell') {
    diffBN_final_result <- run_diffBN(net_struc,data,
                                      meta,index,permutation_time,bootstrap_time,
                                      diffBN_mode)
  } else if (mode == 'bulk') {
    scm.ref = get_scm(mem = data, graph = net_struc, id = "ref") # get reference
    Features <- rownames(data)
    system.time({
      scm.inter = foreach::foreach(i = 1:nrow(data), .combine = cbind) %dopar% {
        get_scm(mem = data[-i,], graph = net_struc, id = Features[i])
      }
    }) # get diffBN for every gene
    diffBN_final_result <- get_diffBN(scm.inter = scm.inter, scm.ref = scm.ref) # ger diffBN summary
  } else {diffBN_final_result=NULL; print('Wrong mode!')}
  return(diffBN_final_result)
}

#' Title get_diffBN
#'
#' get_diffBN: based on a set of bn.fit results and a bn.fir reference, calculate the diffBN results
#'
#' @param scm.inter data.frame of bn.fit results
#' @param scm.ref data.frame of bn.fit reference, which don't have any gene deletion or permutation
#' @param mode character, determine how the diffBN is calculate, usually use "mean"
#'
#' @return data.frame of diffBN
#' @export
#'
#' @examples
#' library(foreach)
#' library(parallel)
#' library(doParallel)
#' registerDoParallel(8) # use 8 cores to speed up the calculation
#' scm.ref = get_scm(mem = ma_test, graph = e, id = "ref") # e is the result of function "trimDAG", is the final netword structure
#' Features <- rownames(ma_test)
#' system.time({
#'   scm.inter = foreach(i = 1:nrow(ma_test), .combine = cbind) %dopar% {
#'     get_scm(ma_test = ma_test[-i,], graph = graph, id = Features[i])
#'   }
#' })
#' diffBN <- get_diffBN(scm.inter = scm.inter, scm.ref = scm.ref)
get_diffBN <- function(scm.inter=NULL, scm.ref=NULL,mode='mean'){
  if (mode=='mean')
  {
    nGene <- ncol(scm.inter)
    diffBN_trans <- rep(scm.ref, nGene) - scm.inter
    diffBN <- diffBN_trans %>% t %>% as.data.frame
  }
  else {
    diffBN=NULL
    warning('Invalid mode!')
  }
  return(diffBN)
}

#' Title run_diffBN
#'
#' run_diffBN is a function that calculate diffBN results based on a single cell expression data
#'
#' @param net_struc bn.fit structure
#' @param data data.frame of single cell expression data
#' @param meta character of metadata corresponding to expression data
#' @param index list of sampling results
#' @param permutation_time integer, times of gene permutation, usually 20-50
#' @param bootstrap_time integer, times of cell sampling, usually 20-50
#' @param diffBN_mode character, determine the diffBN function, usually "mean"
#'
#' @return diffBN results
#' @export
#'
#' @examples
#' library(foreach)
#' library(parallel)
#' library(doParallel)
#' registerDoParallel(8) # use 8 cores to speed up the calculation
#' index <- bootstrap_index(dat_fibro_meta,bootstrap_times=20,ratio=0.4)
#' test_result <- run_diffBN(net_struc,dat_fibro,dat_fibro_meta,index,permutation_time=20,bootstrap_time=20,diffBN_mode='mean')
run_diffBN <- function(net_struc,data,meta,index,permutation_time=20,bootstrap_time=20,diffBN_mode='mean')
{
  diffBN_result <- list(NULL);diff_num=1
  diffBN_raw <- list(NULL)
  diffBN_ref <- list(NULL)
  meta=as.matrix(meta)
  Features <- rownames(data)
  for (i in 1:bootstrap_time)
  {
    data_tmp=data[,index[[i]]]
    if(class(data_tmp)[1]=='Seurat'){data_tmp=as.data.frame(data_tmp@assays$RNA@data)}
    meta_tmp=meta[index[[i]]]
    ds <- gem2mem(data_tmp,meta_tmp,'mean')
    scm.ref = get_scm(mem = ds, graph = net_struc, id = "ref")
    diffBN_ref[[i]] = scm.ref; names(diffBN_ref)[i] = paste0('sample:',i,'_ref')
    for (h in 1:permutation_time)
    {
      scm.inter = foreach::foreach(j = 1:nrow(data_tmp), .combine = cbind) %dopar% {
        get_scm(mem = gem2mem(g_per(data_tmp,j),meta_tmp,'mean'), graph = net_struc, id = Features[j])
      }
      diffBN <- get_diffBN(scm.inter = scm.inter, scm.ref = scm.ref, mode = diffBN_mode)
      diffBN_result[[diff_num]] <- diffBN
      diffBN_raw[[diff_num]] <- scm.inter
      names(diffBN_raw)[diff_num] <- paste0('sample_',i,'/',bootstrap_time,'_permute_',h,'/',permutation_time)
      names(diffBN_result)[diff_num] <- paste0('sample_',i,'/',bootstrap_time,'_permute_',h,'/',permutation_time);diff_num=diff_num+1
      print(paste0('bootstrap_time:',i,'/permutation_time:',h))
    }
  }
  diffBN_final_result <- list(diff=diffBN_result,raw=diffBN_raw,ref=diffBN_ref,permutation_time=permutation_time,bootstrap_time=bootstrap_time)
  return(diffBN_final_result)
}

#' Title
#'
#' @param data_per data.frame
#' @param gene_pos interger
#'
#' @return data.frame
#' @export
#'
#' @examples
#'print('')
g_per <- function(data_per,gene_pos)
{
  data_per[gene_pos,] <- as.numeric(sample(data_per[gene_pos,],length(data_per[gene_pos,])))
  return(data_per)
}

#' Title get_scm
#'
#' get_scm:calculate bn.fit results according to a certain neetwork structure and a cell type expression data.frame
#'
#' @param mem data.frame of cell type expression
#' @param graph bn.fit structure
#' @param id character, id should be given to label the information of sampling and permutation, best in the form of s1p1
#'
#' @return bn.fit results
#' @export
#'
#' @examples
#' scm.ref = get_scm(mem = ma_test, graph = e, id = "ref") # e is the result of function "trimDAG", is the final netword structure
get_scm <- function(mem=NULL, graph=NULL, id=NULL){
  if(is.null(mem)) stop('MEM data is missing')
  if(is.null(graph)) stop('DAG used for linear regression is missing')
  if(is.null(id)) stop('id should be given to label the information of sampling and permutation, best in the form of s1p1')
  
  fit_result <- bnlearn::bn.fit(graph, mem)
  scm <- bnlearn::arcs(graph) %>% as.data.frame
  scm[id] <- sapply(rownames(scm),
                    function(x) stats::coef(fit_result)[[as.character(scm[x,'to'])]][as.character(scm[x,'from'])])
  scm <- scm %>% dplyr::mutate(edge=paste(from,to,sep='~')) %>%
    tibble::column_to_rownames(var = 'edge') %>%
    dplyr::select(-c('from','to'))
  return(scm)
}
