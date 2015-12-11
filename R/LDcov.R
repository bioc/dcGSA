#' Calculate longitudinal distance covariance statistics.
#' @param x.dist A block-diagonal distance matrix of each block being
#'        pairwise distance matrix of genes for each subject.
#' @param y.dist A block-diagonal distance matrix of each block being
#'        pairwise distance matrix of clinical outcomes for each subject.
#' @param nums A vector of integer numbers indicating the number of
#'        repeated measures for each subject.
#' @param bmat A numerical matrix with one column for each subject
#'        (binary values indicating the locations of the repeated measures
#'         for that subject).
#' @return returns the longitudinal distance covariance statistics.
#' @export
#' @examples
#' \dontrun{require(Matrix)}
#' x <- cbind(rnorm(7),rnorm(7)) ## two genes
#' y <- cbind(rnorm(7),rnorm(7)) ## two clinical outcomes
#' ## Two subjects: the first one has three measures
#' ## while the other one has four measures
#' ID <- c(1,1,1,2,2,2,2) ## The IDs for the two subjects.
#' nums <- c(3,4) ## number of repeated measures for each subjects
#' ## prepare block-diagonal distance matrix for genes and clinical outcomes
#' lmat <- lapply(nums,function(x){z=matrix(1,nrow=x,ncol=x)})
#' mat <- as.matrix(bdiag(lmat))
#' lmat <- lapply(nums,function(x){z=matrix(0,nrow=x,ncol=x);z[,1]=1;z})
#' bmat <- as.matrix(bdiag(lmat))
#' ind <- apply(bmat,2,sum)
#' bmat <- bmat[,ind!=0]
#' ydist <- as.matrix(dist(y))*mat
#' xdist <- as.matrix(dist(x))*mat
#
#' LDcov(x.dist=xdist,y.dist=ydist,nums=nums,bmat)

LDcov <- function(x.dist=NULL,y.dist=NULL,nums=NULL,bmat=NULL){
  s1 <- apply((x.dist*y.dist)%*%bmat,2,sum)/nums^2
  s2 <- apply((x.dist%*%bmat),2,sum)*apply((y.dist%*%bmat),2,sum)/nums^4
  n=ncol(x.dist)
  s3 <- apply(cbind(x.dist,y.dist),1,function(x){
    sum(x[1:n]%*%t(x[(n+1):(2*n)]))
    })
  s3 <- as.numeric(s3%*%bmat)/nums^3
  mean(nums*(s1+s2-2*s3)/s2,na.rm=TRUE)
}
