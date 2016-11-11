#' Read gene set file in gmt format
#' @param file filename for the gmt file
#' @return a list of gene sets with each element being a vector of gene names
#' @export
#' @examples
#' fpath <- system.file("extdata", "sample.gmt.txt", package="dcGSA")
#' GS <- readGMT(file=fpath)
readGMT <- function(file = NULL){
  geneset <- read.table(file=file,header=FALSE,sep="\n")
  geneset <- as.character(geneset[,1])
  geneset<- sapply(geneset,function(x){strsplit(x,split="\t")})
  names(geneset) <- sapply(geneset,function(x){x[1]})
  geneset <- lapply(geneset,function(x){x[-(1:2)]})
  return(geneset)
}
