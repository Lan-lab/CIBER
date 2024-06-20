#' Generate edge set
#'
#' @param fromSet list of starting nodes
#' @param toSet list of ending nodes
#' @param sep character, used for sepearating names
#'
#' @return list of edge names
#' @export
#'
#' @examples
#' ctypes <- c("MPPa", "MPPb")
#' s.diffScore <- setEdges(ctypes, ctypes, sep = "~")
setEdges <- function(fromSet, toSet, sep = "~") {
  edges <- paste(rep(fromSet, each = length(toSet)), toSet, sep = sep)
  selfs <- paste(fromSet, fromSet, sep = sep)
  return(setdiff(edges, selfs))
}

#' Get diffScore
#'
#' @param diffBN data frame, with diffBN data previously calculated
#' @param edgeSet list of edge names
#' @param abs logical, set to TRUE when diffScore is considered the sum of the absolute values of diffBN
#'
#' @return list of diffScores
#' @export
#'
#' @examples
#' ds.diffScore <- diffScore(diffBN, s.diffScore, abs = TRUE)
diffScore <- function(diffBN, edgeSet, abs = TRUE) {
  if (abs) {
    diffBN <- abs(diffBN)
  }
  dm <- diffBN[intersect(colnames(diffBN), edgeSet)]
  return(rowSums(dm))
}

#' Rank genes according to diffScores
#'
#' @param diffScores list of previously calculated diffScores
#'
#' @return list of ranked gene names
#' @export
#'
#' @examples
#' gl.diffScore <- dsRank(ds.diffScore)
dsRank <- function(diffScores) {
  diffScores %>%
    sort(decreasing = TRUE) %>%
    names()
}
