#' Read in binary files containing genomic scoring results
#' 
#' @param genome TO_BE_ADDED
#' @param score.dir TO_BE_ADDED
#' @param prot.name TO_BE_ADDED
#' @param chr.name TO_BE_ADDED
#' @param k TO_BE_ADDED
#'
#' @return TO_BE_ADDED
#' 
#' @examples
#'
#' @export
#' 
score.parser = function(genome, score.dir, prot.name, chr.name, k) {
    s = length(genome[[chr.name]])
    fwd.score = readBin(con = paste0(score.dir,prot.name,"_",chr.name,"_F.dat"), what = "double", n=s, endian = "big")
    rev.score = readBin(con = paste0(score.dir,prot.name,"_",chr.name,"_R.dat"), what = "double", n=s, endian = "big")
    flag = readBin(con = paste0(score.dir,prot.name,"_",chr.name,"_Flag.dat"), what = "integer", n=s, size=1, endian = "big")
    return(list(score=rbind(fwd.score, rev.score), flag=flag, chr=chr.name, k=k))
  }
