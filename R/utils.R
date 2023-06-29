#' Drop rows with missing values or zero values
#'
#' This function drops rows in a matrix that have missing values, but keeps values equal to zero.
#'
#' @param x A matrix
#' @return A matrix with NA dropped
#' @export
dropNA <- function(x) {
  if(!is(x, "matrix")) stop("x needs to be a matrix!")
  
  zeros <- which(x==0, arr.ind=TRUE)
  ## keep zeros
  x[is.na(x)] <- 0
  x[zeros] <- NA
  x <- Matrix::drop0(x)
  x[zeros] <- 0
  x
}

#' Minmax normalization
#'
#' This function normalizes a numeric vector to the range [0,1] using min-max normalization.
#'
#' @param x a numeric vector to be normalized.
#' @return a numeric vector with values between 0 and 1.
#'
minmax.normalisation <- function(x)
{
  return((x-min(x,na.rm=T)) /(max(x,na.rm=T)-min(x,na.rm=T)))
}

#' Sort absolute
#'
#' This function sorts a data table based on the absolute values of a specified column.
#'
#' @param dt a data table to be sorted.
#' @param sort.field a character string indicating the name of the column to sort based on its absolute value.
#' @return a sorted data table.
#'
sort.abs <- function(dt, sort.field) dt[order(-abs(dt[[sort.field]]))]

#' Turn data table to matrix
#'
#' This function converts a data table (excluding the first column) to a matrix and sets the row names to the values of the first column.
#'
#' @param x a data table to be converted to a matrix.
#' @return a matrix with row names taken from the first column of the input data table.
#'
matrix.please<-function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}

#' Annotates peaks with genes in gene window
#'
#' @param gene_metadata.dt a data table of gene metadata
#' @param atac.sce a single-cell experiment object
#' @param distance a numeric value specifying the distance from the gene body to the peak
#' 
#' @import data.table
#' @import dplyr
#' @importFrom stringr str_split
#' @importFrom purrr map_chr
#' 
#' @return a data table with the overlap between peaks and genes, including the distance between the peak and the gene body
#' 
#' @export
annotate_peaks = function(gene_metadata.dt = gene_metadata.dt,
                          atac.sce = atac.sce,
                          distance = 1.5e5){
    
    # Load gene metadata
    gene_metadata <- copy(gene_metadata.dt) %>% 
      .[,chr:=as.factor(sub("chr","",chr))] %>%
      setnames("symbol","gene") %>%
      .[, c("chr","start","end","gene","ens_id","strand")]

    peakSet.dt <- data.table(peak = rownames(atac.sce)) %>%
        .[,`:=`(chr = as.factor(peak %>% str_split(':') %>% map_chr(1) %>% gsub('chr', '', .)),
                start = as.integer(peak %>% str_split(':') %>% map_chr(2) %>% str_split('-') %>% map_chr(1)),
                end = as.integer(peak %>% str_split(':') %>% map_chr(2) %>% str_split('-') %>% map_chr(2)))] %>%
        .[,c('chr', 'start', 'end', 'peak')] %>% 
      setkey(chr,start,end)

    ## Overlap 
    gene_metadata.ov <- copy(gene_metadata) %>%
      .[strand=="+",c("gene.start","gene.end"):=list(start,end)] %>%
      .[strand=="-",c("gene.start","gene.end"):=list(end,start)] %>%
      .[strand=="+",c("start","end"):=list (gene.start-gene_window, gene.end+gene_window)] %>%
      .[strand=="-",c("end","start"):=list (gene.start+gene_window, gene.end-gene_window)] %>% 
      setkey(chr,start,end)

    stopifnot((gene_metadata.ov$end-gene_metadata.ov$start)>0)

    ov <- foverlaps(
      peakSet.dt,
      gene_metadata.ov,
      nomatch = NA
    ) %>%  .[,c("start","end"):=NULL] %>%
      setnames(c("i.start","i.end"),c("peak.start","peak.end")) %>%
      .[,peak.mean:=(peak.start+peak.end)/2] %>%
      # calculate distance from the peak to the genebody
      .[,dist:=min(abs(gene.end-peak.mean), abs(gene.start-peak.mean)), by=c("gene","ens_id","peak","strand")] %>%
      .[strand=="+" & peak.mean>gene.start & peak.mean<gene.end,dist:=0] %>%
      .[strand=="-" & peak.mean<gene.start & peak.mean>gene.end,dist:=0]
    
    return(ov)   
}
