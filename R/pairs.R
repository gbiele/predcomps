#' GetPairs
#'
#' Form all pairs of rows in \code{X} and compute Mahalanobis distances based on \code{v}. 
#' 
#' To help with computational constraints, you have the option to not form pairs between all rows of \code{X} but instead of specify a certain number (\code{numForTransitionStart}) to randomly be selected as rows from which transitions start, and another number (\code{numForTransitionEnd}) to be randomly selected as where transitions end. We then form all pairs between transition-start rows and transition-end rows. 
#' 
#' In order to get a smaller data frame for later manipulations (and maybe just because it's a good idea), you can also specify \code{onlyIncludeNearestN}, in which case we return only the nearest \code{onlyIncludeNearestN} transition ends for each transition start (instead of all pairs).
#' 
#' @param X data frame
#' @param u input of interest
#' @param v other inputs
#' @param mahalanobisConstantTerm Weights are (1 / (mahalanobisConstantTerm + Mahalanobis distance))
#' @param numForTransitionStart number of rows to use as the start points of transitions (defaulting to `NULL`, we use all rows)
#' @param numForTransitionEnd number of rows to use as potential end points of transitions (defaulting to `NULL`, we use all rows)
#' @param onlyIncludeNearestN for each transition start, we only include as transition end points the nearest `onlyIncludeNearestN` rows  (defaulting to `NULL`, we use all rows)
#' @return a data frame with the inputs \code{v} from the first of each pair, \code{u} from each half (with ".B" appended to the second), and the Mahalanobis distances between the pairs.
#' @export
#' @examples
#' v <- rnorm(100)
#' u <- v + 0.3*rnorm(100)
#' qplot(v,u)
#' X = data.frame(v=v,u=u)
#' pairsDF <- GetPairs(X, "v", "u")
#' pairsDFRow1 <- subset(pairsDF, OriginalRowNumber==1)
#' # When we subset to one "original row number", all of the v's are the same:
#' print(pairsDFRow1$v)
#' # ... and u's corresponding to closer v.B (the v in the second element of the pair) have higher weight:
#' qplot(u.B, Weight, data=pairsDFRow1)
GetPairs <- function(X, u, v,
                     numForTransitionStart = NULL,
                     numForTransitionEnd = NULL,
                     onlyIncludeNearestN = NULL,
                     mahalanobisConstantTerm=1) {
  
  assert_that(length(u) == 1) # make sure we have exactly 1 input var of interest
  for (columnName in c(u,v)) {
    assert_that(columnName %in% names(X))
    columnClass <- class(X[[columnName]])
    if (!(columnClass) %in% c("integer", "numeric")) {
      stop(sprintf("Sorry, column %s is of class %s. I can only deal with integer and numeric types for now.", columnName, columnClass))
    }
  }
  ######################## start updated section (GB) ###############################
  if (!is.null(numForTransitionStart)) {
  # I just named rows like this to be sure that I calculate distances
  # for the sampling case correctly
    s1 <- sample.int(nrow(X), size=numForTransitionStart) 
    X1 <- X[s1, c(v,u), with = F]
    X1$OriginalRowNumber <- s1
  } else {
    X1 <- X[,c(v,u), with = F]
    s1 <- 1:1:nrow(X1)
  }
  X1$OriginalRowNumber <- s1
  
  if (!is.null(numForTransitionEnd)) {
    s2 <- sample.int(nrow(X), size=numForTransitionEnd)
    X2 <- X[s2, c(v,u), with = F]
    X2$OriginalRowNumber.B <- s2
  } else {
    X2 <- X[,c(v,u), with = F]
    s2 <- 1:nrow(X)
  }
  X1$OriginalRowNumber <- s1
  X2$OriginalRowNumber.B <- s2
  
  
  vMatrix1 <- as.matrix(X1[,v, with = F])
  vMatrix2 <- as.matrix(X2[,v, with = F])
  vMatrix = as.matrix(X[unique(c(s1,s2)),v,with = F])
  covV=cov(vMatrix)
  # because memory is an issue for my datasets, 
  # I delete variabeles if they are no longer needed
  rm(vMatrix) 
  
  # the maha function from mvnfast is much faster than standard R
  # but the majof speed gain is not made here
  distMatrix <- apply(vMatrix1, 1, function(row) maha(vMatrix2, row, covV))
  #dim(distMatrix)
  rm(vMatrix1,vMatrix2,covV)
  colnames(distMatrix) <- s1
  rownames(distMatrix) <- s2
  if (is.null(numForTransitionEnd)) {
    # Here we are saving lots of RAM!
    # instead of creating distDF for the entire data matrix, 
    # distDF is formed only from the relevant lower triangular matrix
    # combnPrim is a fast version of combn, which is relevant for large data sets
    # minor issues: combnPrim comes with gRbase, which can only be installed via bioconductor
    distDF = data.table(t(rbind((combnPrim(1:dim(X1)[1],2)),distMatrix[lower.tri(distMatrix)])))
    } else {
    # now the more complicated solution for sampled pairs
    # this is base R, I haven't checked if this could be speeded up
      distDF = data.table(t(rbind(as.vector(t(matrix(rep(s1,length(s2)),nrow = length(s1),ncol = length(s2)))),
                                  as.vector(matrix(rep(s2,length(s1)),nrow = length(s2),ncol = length(s1))),
                                  as.vector(distMatrix))))
  }
  
  rm(distMatrix)
  setnames(distDF,names(distDF),c("OriginalRowNumber", "OriginalRowNumber.B", "MahalanobisDistance"))
  
  ######################## end updated section (GB) ###############################
  ### there are still some efficientcy gains below because we use data.tables ####
  
  if (!is.null(onlyIncludeNearestN)) {
    distDF <- distDF %.% 
      group_by(OriginalRowNumber) %.% 
      filter(rank(MahalanobisDistance, ties.method="random") < onlyIncludeNearestN)
  }
  # note that beceause X1 and distDF are data tables, 
  # we are calling date.tables merge function here, which is much faster than merging data.frames
  pairs <- merge(X1, distDF, by = "OriginalRowNumber") 
  pairs <- merge(X2, pairs, by = "OriginalRowNumber.B", suffixes = c(".B", ""))
  pairs$Weight <- 1/(mahalanobisConstantTerm + pairs$MahalanobisDistance)
  
  # If we haven't sampled, then OriginalRowNumber == OriginalRowNumber.B means that 
  # the transition start and end are the same, so we should remove those rows.
  if (is.null(numForTransitionStart) && is.null(numForTransitionEnd)) {
    pairs <- subset(pairs, OriginalRowNumber != OriginalRowNumber.B)
  }
  
  # Renormalize weights:
  # this could also maybe be speeded up by using data.table isteadn of dplyr
  # then again, the group_by probably collabses redundant rows (sseing this only now. 13. Osctober)
  # so something like
  # pairs[, Weight := Weight/sum(Weight)]
  pairs <- pairs %.% group_by(OriginalRowNumber) %.% mutate(Weight = Weight/sum(Weight))
  
  return(data.frame(pairs))
}


#' PlotPairCumulativeWeights
#'
#' For a sample of transition start rows, we plot rank of transition end (by increasing weight) vs. cumulative weight. This gives a sense of how much weight is going into the nearest points vs. further ones.
#' 
#' @export
#' @examples
#' v <- rnorm(100)
#' u <- v + 0.3*rnorm(100)
#' X = data.frame(v=v,u=u)
#' pairsDF <- GetPairs(X, "v", "u")
#' pairsDFRow1 <- subset(pairsDF, OriginalRowNumber==1)
#' # For most original rows, we get 75% of the weight in 50% of the pairs:
#' PlotPairCumulativeWeights(pairsDF)

PlotPairCumulativeWeights <- function(pairs, numOriginalRowNumbersToPlot = 20) {
  rowNumSample <- sample(unique(pairs$OriginalRowNumber))[1:numOriginalRowNumbersToPlot]
  pairsWithCumWeightSums <- pairs %.% 
    group_by(OriginalRowNumber) %.% 
    arrange(OriginalRowNumber, -Weight) %.% 
    mutate(CumulativeWeight = cumsum(Weight), Rank = dense_rank(-Weight))
  
  pairsSubset <- subset(pairsWithCumWeightSums, OriginalRowNumber %in% rowNumSample)
  
  ggplot() + 
    geom_line(aes(x=Rank, y=CumulativeWeight, color=factor(OriginalRowNumber)), data = pairsSubset, alpha = .2) + 
    geom_line(aes(x=Rank, y=CumulativeWeight), stat = "summary", fun.y = "median", data=pairsWithCumWeightSums)
}  
