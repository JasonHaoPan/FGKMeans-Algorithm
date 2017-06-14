library(dplyr)
library(wskm)
# Initialize chromosome
initialize_chrom <- function(){
  # We can calculate the length of chromosome from feature number and and featrue group number
  # by multiply
  chrom_len <- feature_num*feature_group_num
  chrom_matr <- matrix(0, num_chrom, chrom_len) 
  
  i <- 1 
  for (i in 1:num_chrom){
    temp_arr <- matrix(0, feature_group_num, feature_num) 
    j <- 1
    while (j <= feature_num){
      #for each chromosome，it can only belong to one feture group
      r=sample(1:feature_group_num, 1) 
      temp_arr[r,j] = 1
      j <- j + 1
    }
    #
    temp_arr <- check_matr(temp_arr) 
    
    chrom <- c(temp_arr)
    k <- 1
    for (k in 1:chrom_len){
      chrom_matr[i,k] <- chrom[k]
      k <- k + 1
    }
  }
  return (chrom_matr)
  
}


check_matr <- function(x){
  for (i in 1:feature_group_num){
    if (sum(x[i,]) == 0){ 
      for (j in 1:feature_group_num){
        if (sum(x[j,]) > 1){ 
          for (k in 1:feature_num){
            if (x[j,k] == 1){ 
              x[j,k] <- 0 
              x[i,k] <- 1
              break 
            }
          }
        }
      }
    }
  }
  return(x)
}

arr_to_group <- function(x){
  i <- 1
  feat_to_group <- c() 
  for (i in 1:length(x)){
    if (x[i] == 1){
      if (i %% feature_group_num == 0) 
        i <- feature_group_num
      else
        i <- i %% feature_group_num
      feat_to_group <- c(feat_to_group,i) 
      groupInfo <<- feat_to_group
      
    }
  }
  return (groupInfo)
}


fgkm <- function(x, centers, groups, lambda, eta, maxiter=100, delta=0.000001, maxrestart=10,seed=-1) 
{
  if (missing(centers))
    stop("the number or initial clusters 'centers' must be provided")
  
  if(seed<=0){
    seed <-runif(1,0,10000000)[1]
  }
  
  vars <- colnames(x)
  
  nr <-nrow(x) # nrow() return a integer type 
  nc <-ncol(x) # integer 
  
  if (is.data.frame(centers) || is.matrix(centers))
  {
    init <- TRUE
    k <- nrow(centers)
  }
  else
  {
    init <- FALSE
    k <- centers
    centers <- double(k * nc)
  }
  
  # get the setting of feature group
  if (is.character(groups) && length(groups) == 1) {
    G <- .C("parseGroup",as.character(groups),numGroups=integer(1), groupInfo=integer(nc),PACKAGE="wskm")
    
  } else if (is.vector(groups) && length(groups) == nc) {
    G <- list()
    grps <- as.factor(groups)
    groupNames <- levels(grps)
    G$numGroups <- nlevels(grps)
    G$groupInfo <- as.integer(as.integer(grps) - 1)
  }
  
  set.seed(seed)
  Z <- .C("fgkm",
          x = as.double(as.matrix(x)),
          nr,
          nc,
          k = as.integer(k),
          lambda = as.double(lambda),
          eta = as.double(eta),
          G$numGroups,
          G$groupInfo,
          delta = as.double(delta),
          maxIterations = as.integer(maxiter),
          maxRestarts = as.integer(maxrestart),
          as.logical(init),
          #          seed,
          cluster = integer(nr),
          centers=as.double(as.matrix(centers)),
          featureWeight = double(k * nc),
          groupWeight = double(k * G$numGroups),
          iterations = integer(1),
          restarts = integer(1),
          totiters = integer(1),
          totalCost = double(1),
          totss = double(1),
          withinss = double(k),
          PACKAGE="wskm"
  )
  
  centers <- matrix(Z$centers)
  dim(centers) <- c(k, nc)
  colnames(centers) <- vars
  
  featureWeight <- matrix(Z$featureWeight)
  dim(featureWeight) <- c(k, nc)
  colnames(featureWeight) <- vars
  
  groupWeight <- matrix(Z$groupWeight)
  dim(groupWeight) <- c(k, G$numGroups )
  colnames(groupWeight) <- 1:ncol(groupWeight)
  
  ignore <- which(rowSums(centers==0) == ncol(centers))
  if (length(ignore)) {
    centers <- centers[-ignore,, drop=FALSE]
    featureWeight <- featureWeight[-ignore,, drop=FALSE]
  }
  
  rownames(centers) <- 1:nrow(centers)
  rownames(featureWeight) <- 1:nrow(featureWeight)
  rownames(groupWeight) <- 1:nrow(groupWeight)
  
  cluster <- Z$cluster + 1
  
  size <- aggregate(cluster, list(cluster=cluster), length)[[2]]
  
  result <- list(cluster = cluster,
                 centers = Z$centers,
                 totss = Z$totss, 
                 withinss = Z$withinss, 
                 tot.withinss = sum(Z$withinss), 
                 betweenss = Z$totss-sum(Z$withinss),
                 size = size,
                 iterations = Z$iterations,
                 restarts = Z$restarts,
                 totiters=Z$totiters,
                 featureWeight = Z$featureWeight,
                 groupWeight = Z$groupWeight)
  
  dim(result$centers) <- c(k, nc)
  dim(result$featureWeight) <- c(k, nc)
  dim(result$groupWeight) <- c(k, G$numGroups)
  
  class(result) <- c("kmeans", "fgkm")
  return(result)
}

test_data <- read.csv("~/Desktop/SZU/1_Dataset/Lymphoma.csv", header = TRUE, sep = ",")
# print(test_data)
test_data <- test_data[,-1]
lambda= 1
eta = 1

object_num <<- nrow(test_data)
#print object number
print(object_num)

feature_num <<- ncol(test_data) 
#print feature number
print(feature_num)
feature_group_num = 10
cluster_num = 3
num_chrom = 20 
chroms <- initialize_chrom()
#print(chroms)
centers = 3
lst = list()
for(i in 1:num_chrom){
  groups <- arr_to_group(chroms[i,])
  temp = fgkm(test_data, centers, groups, lambda, eta)
  lst[i] = list(temp)
}



