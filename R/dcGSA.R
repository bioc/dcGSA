#' Perform gene set analysis for longitudinal gene expression profiles.
#' @param data A list with ID (a character vector for subject ID),
#'        pheno (a data frame with each column being one clinical outcome),
#'        gene (a data frame with each column being one gene).
#' @param geneset A list of gene sets of interests (the output of
#'        \code{\link{readGMT}} function).
#' @param nperm An integer number of permutations performed to get P values.
#' @param c An integer cutoff value for the overlapping number of genes
#'        between the data and the gene set.
#' @param KeepPerm A logical value indicating if the permutation statistics are
#'         kept. If there are a large number of gene sets and the number of 
#'         permutation is large, the matrix of the permutation statistics could
#'         be large and memory demanding.
#' @param parallel A logical value indicating if parallel computing is wanted.
#' @param BPparam Parameters to configure parallel evaluation environments
#'         if parallel is TRUE. The default value is to use 4 cores in a single
#'         machine. See \code{\link{BiocParallelParam}} object in Bioconductor
#'         package \code{BiocParallel} for more details.
#' @return Returns a data frame with following columns, if KeepPerm=FALSE;
#'         otherwise, returns a list with two objects: "res" object being the 
#'         following data frame and "stat" being  the permutation statistics.
#' \item{Geneset}{Names for the gene sets.}
#' \item{TotalSize}{The original size of each gene set.}
#' \item{OverlapSize}{The overlapping number of genes between
#'                    the data and the gene set.}
#' \item{Stats}{Longitudinal distance covariance between the clinical outcomes
#'              and the gene set.}
#' \item{NormScore}{Only available when permutation is performed.
#'       Normalized longitudinal distance covariance using the mean and
#'       standard deviation of permutated values.}
#' \item{P.perm}{Only available when permutation is performed. Permutation P 
#'               values.}
#' \item{P.approx}{P values obtained using normal distribution to approximate 
#'                 the null distribution.}
#' \item{FDR.approx}{FDR based on the P.approx.}
#' @references Distance-correlation based Gene Set Analysis in Longitudinal
#'             Studies. Jiehuan Sun, Jose Herazo-Maya, Xiu Huang,
#'             Naftali Kaminski, and Hongyu Zhao.
#' @export
#' @examples
#' data(dcGSAtest)
#' fpath <- system.file("extdata", "sample.gmt.txt", package="dcGSA")
#' GS <- readGMT(file=fpath)
#' system.time(res <- dcGSA(data=dcGSAtest,geneset=GS,nperm=100))
#' head(res)
dcGSA <- function(data=NULL,geneset=NULL,nperm=10,c=0,KeepPerm=FALSE,
                  parallel=FALSE,BPparam = MulticoreParam(workers=4)){

  ord <- order(data$ID)
  data$ID <- as.character(data$ID[ord])
  data$pheno <- as.data.frame(data$pheno)[ord,]
  data$pheno <- as.data.frame(data$pheno)
  data$gene <- as.data.frame(data$gene)[ord,]
  data$gene <- as.data.frame(data$gene)

  ids <- unique(as.character(data$ID))
  id <- factor(as.character(data$ID),levels=ids)
  nums <- sapply(ids,function(x){sum(data$ID==x)})
  totalsize <- sapply(geneset, function(x){length(x)})
  geneset.list <- lapply(geneset,function(x)
    {x[!is.na(match(x,colnames(data$gene)))]
    })
  size <- sapply(geneset.list, function(x){length(x)})

  if(sum(size > c) < 1){
    stop(paste("No geneset has overlapping size greater than ",c,"!",sep="") )}
  if(sum( nums < 3 ) > 0){
    stop("each subject must have more than two visits!")}
  if(sum( sapply(data,function(x) sum(is.na(x))) ) > 0 ){
    stop("no missing values are allowed!")}

  lmat <- lapply(nums,function(x){z=matrix(1,nrow=x,ncol=x)})
  mat <- as.matrix(bdiag(lmat))
  bmat <- model.matrix( ~ -1+id,
                        contrasts.arg  = list(id=contrasts(id,contrasts=FALSE)))

  geneset.na <- names(geneset.list)[size <= c]
  geneset.list <- geneset.list[size > c]
  totalsize <- totalsize[size > c]
  size <- size[size > c]

  ydist <- as.matrix(dist(data$pheno))*mat

  xdist.list <- lapply( geneset.list, function(x){
    as.matrix(dist(data$gene[x]))*mat
  })

  stats <- sapply(xdist.list,LDcov,y.dist=ydist,nums=nums,bmat=bmat)
  res <- data.frame(Geneset=names(geneset.list),TotalSize=totalsize,
                    OverlapSize=size,Stats=stats)

  if(nperm > 0){

    y.split <- split(data$pheno,as.character(data$ID))
    if(parallel){
      res.per <- bplapply(seq_len(nperm),function(i){
        y.per<- lapply(seq_along(nums),function(x){
          as.data.frame(y.split[[ids[x]]][sample(1:nums[x],nums[x]),])
          })
        ydist.per <- as.matrix(dist(do.call(rbind,y.per)))*mat
        sapply(xdist.list,LDcov,y.dist=ydist.per,nums=nums,bmat=bmat)
      },BPPARAM=BPparam)
      res.per <- do.call(cbind,res.per)
    }else{
      res.per <- lapply(seq_len(nperm),function(i){
        y.per<- lapply(seq_along(nums),function(x){
          as.data.frame(y.split[[ids[x]]][sample(seq_len(nums[x]),nums[x]),])
          })
        ydist.per <- as.matrix(dist(do.call(rbind,y.per)))*mat
        sapply(xdist.list,LDcov,y.dist=ydist.per,nums=nums,bmat=bmat)
      })
      res.per <- do.call(cbind,res.per)
    }

    p.val <- sapply(seq_len(nrow(res.per)), function(i){
      (sum(abs(res.per[i,]) >= abs(res$Stats[i])) + 1)/(nperm+1)
    })
    sd.per <- apply(res.per,1,sd)
    mean.per <- apply(res.per,1,mean)
    norm.score <- (res$Stats-mean.per)/sd.per

    res <- data.frame(Geneset=names(geneset.list),TotalSize=totalsize,
                      OverlapSize=size,Stats=stats,NormScore=norm.score,
                      P.perm = p.val, P.approx = 1-pnorm(norm.score),
                      FDR.approx = p.adjust(1-pnorm(norm.score),"fdr"))
    res <- res[order(res$NormScore,1-res$P.perm,decreasing = TRUE),]
    rownames(res) <- NULL
    if(KeepPerm){
        res = list(res=res,stat=res.per)
    }
    res
  }else{
    res <- res[order(res$Stats,decreasing = TRUE),]
    rownames(res) <- NULL
    res
  }
}
