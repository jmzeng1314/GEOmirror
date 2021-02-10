##' Download from geoChina by a GSE id
##'
##' \code{geoChina} will download the expression matrix and phenotype data as an object in ExpressionSet format
##'      from cloud in mainland China, so it will be a little fast for users in mainland China.
##'      So far, only tencent mirror for our geoChina.
##'      it's just a alternative method for getGEO function from GEOquery package.
##'      geoChina('GSE1009') is the same as eSet=getGEO('GSE1009', getGPL = F)
##'
##' @param gse input GSE id, such as GSE1009, GSE2546, default:GSE2546
##' @import GEOquery
##' @return a list of ExpressionSet, which contains the  expression matrix and phenotype data
##' @examples
##' geoChina()
##' geoChina('GSE1009')
##' geoChina('GSE2546')
##' @export geoChina

geoChina <- function(gse='GSE2546',mirror='tencent'){
  # same as code : eSet=getGEO('GSE2546', destdir=".", AnnotGPL = F, getGPL = F)
  # http://49.235.27.111/GEOmirror/GSE2nnn/GSE2546_eSet.Rdata
  # gse='GSE2546';mirror='tencent'

  gse=toupper(gse)

  if(!gse %in% series.accession){
    stop('Your GSE may not be expression by array, or even not a GSE')
  }

  down=ifelse(as.numeric(gsub('GSE','',gse))<1000,
              paste0('/GEOmirror/GSEnnn/',gse,
                     '_eSet.Rdata'),
              paste0('/GEOmirror/',
                     gsub('[0-9][0-9][0-9]$','nnn',gse),'/',gse,
                     '_eSet.Rdata'))

  if(mirror=='tencent'){
    up='http://49.235.27.111'
  }

  tpf=paste0(gse, '_eSet.Rdata')
  download.file(paste0(up,down),tpf,mode = "wb")
  suppressWarnings(load(tpf))

  # getGEO('GSE2546', destdir=".", AnnotGPL = F, getGPL = F)
  cat(paste0('you can also use getGEO from GEOquery, by \ngetGEO(',
             shQuote(gse),
             ', destdir=".", AnnotGPL = F, getGPL = F)'
  ))
  return(gset)

}
