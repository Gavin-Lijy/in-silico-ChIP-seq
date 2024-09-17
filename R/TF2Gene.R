#' TF2Gene function: Linking TFs to target genes
#'
#' This function takes as input an RNA-Seq dataset, a motif-matching scores dataset, and a dataset linking
#' genomic peaks with genes. It performs a regression analysis to determine the influence of each transcription
#' factor on each gene. The output is a data frame of regression coefficients (beta), p-values, and TFs and genes.
#'
#' @param rna.sce A SingleCellExperiment object containing RNA-seq data.
#' @param motifmatcher_chip.se A MotifMatchR object (output of silico_chip function())
#' @param peak2gene.dt Data table linking genomic peaks with genes
#' @param gene_metadata.dt (optional, Only needed to find peak-genes if peak2gene.dt not supplied) gene metadata used to annotate peaks with genes
#' @param atac.sce (optional, Only needed to find peak-genes if peak2gene.dt not supplied) A SingleCellExperiment object containing the ATAC-seq data
#' @param filter_genes logical indicating whether to filter out non-informative genes (default: TRUE)
#' @param min_chip_score minimum ChIP score for a peak to be considered a target of a TF (default: 0.15)
#' @param distance maximum distance (in bp) between a peak and a gene for them to be considered linked (default: 5e4, can not be larger than maximum distance in peak2gene.dt)
#' @param cores number of cores to use for parallel processing (default: detectCores())
#'
#' @import data.table
#' @import dplyr
#' @import SingleCellExperiment
#' @importFrom parallel detectCores mclapply
#' @importFrom stats lm 
#' @importFrom stringr str_to_title
#'
#' @return a data frame of regression coefficients (beta), p-values, and TFs and genes
#'
#' @export

TF2Gene = function(rna.sce = rna.sce, 
                   motifmatcher_chip.se = motifmatcher_chip.se,
                   peak2gene.dt = peak2gene.dt,
                   gene_metadata.dt = NULL, # Only needed to find peak-genes if peak2gene.dt not supplied
                   atac.sce = NULL, # Only needed to find peak-genes if peak2gene.dt not supplied
                   filter_genes = TRUE,
                   min_chip_score = 0.15,
                   distance = 5e4,
                   cores = detectCores()){
    
    if(is.null(peak2gene.dt)){
        if(is.null(gene_metadata.dt) | is.null(atac.sce)){
            cat('Please provide gene metadata and atac.sce \n')
            stop()
        }
        cat('Annotating peaks with genes in distance window \n')
        peak2gene.dt = annotate_peaks(gene_metadata.dt = gene_metadata.dt,
                                  atac.sce = atac.sce,
                                  distance = distance)
    }
   peak2gene.dt = peak2gene.dt  %>%
                  .[dist<=distance] %>%
                  .[,peak:=sprintf("chr%s:%s-%s",chr,peak.start,peak.end)]
    
    virtual_chip.mtx = assay(motifmatcher_chip.se, 'VirtualChipScores')

    # Sanity checks
    stopifnot(length(intersect(rownames(virtual_chip.mtx),unique(peak2gene.dt$peak)))>1e5)

    ## Link TFs to target genes using the virtual ChIP-seq 
    tf2gene_chip.dt <- mclapply(colnames(virtual_chip.mtx), function(i){
    # Select target peaks (note that we only take positive correlations into account)
    target_peaks_i <- names(which(virtual_chip.mtx[,i]>=min_chip_score))

    if (length(target_peaks_i)>=1) {
        tmp <- data.table(
          tf = i,
          peak = target_peaks_i,
          chip_score = virtual_chip.mtx[target_peaks_i,i]
        ) %>% merge(peak2gene.dt[peak %in% target_peaks_i,c("peak","gene","dist")], by="peak")
        return(tmp)
        }
    }, mc.cores=cores) %>% data.table::rbindlist(.)     
                       
    ## Filter TFs and genes 

    TFs <- intersect(unique(tf2gene_chip.dt$tf),toupper(rownames(rna.sce)))
    genes <- intersect(unique(tf2gene_chip.dt$gene),rownames(rna.sce))
    
    if(filter_genes){
        # filter out non-informative genes
        genes <- genes[grep("*Rik|^Gm|^mt-|^Rps|^Rpl", genes, invert=T)] 
    }

    tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TFs & gene%in%genes,]

    # Fetch RNA expression matrices
    rna_tf.mtx <- assay(rna.sce,"logcounts")[unique(tf2gene_chip.dt$tf),,drop=F]; rownames(rna_tf.mtx) <- toupper(rownames(rna_tf.mtx))
    rna_targets.mtx <- assay(rna.sce,"logcounts")[unique(tf2gene_chip.dt$gene),,drop=F]

    # Filter out lowly variable genes and TFs
    rna_tf.mtx <- rna_tf.mtx[apply(rna_tf.mtx,1,var)>=0.1,,drop=F]
    rna_targets.mtx <- rna_targets.mtx[apply(rna_targets.mtx,1,var)>=0.1,]

    TFs <- intersect(unique(tf2gene_chip.dt$tf),rownames(rna_tf.mtx))
    genes <- intersect(unique(tf2gene_chip.dt$gene),rownames(rna_targets.mtx))
    tf2gene_chip.dt <- tf2gene_chip.dt[tf%in%TFs & gene%in%genes,]

    cat(sprintf("Number of TFs: %s \n",length(TFs)))
    cat(sprintf("Number of genes: %s \n",length(genes)))
    
    ## run regression 
    GRN_coef.dt = mclapply(genes, function(i){
        tfs <- unique(tf2gene_chip.dt[gene==i,tf])
        tfs %>% map(function(j) {
          x <- rna_tf.mtx[j,]
          y <- rna_targets.mtx[i,]
          lm.fit <- lm(y~x)
          data.frame(tf=j, gene=i, beta=round(coef(lm.fit)[[2]],3), pvalue=format(summary(lm.fit)$coefficients[2,4], digits=3))
        }) %>% data.table::rbindlist(.)
      }, mc.cores=cores) %>% data.table::rbindlist(.)
    
    # Add chip score back in so all info is 
    #chip_GRN = merge(tf2gene_chip.dt, GRN_coef.dt, by=c('tf', 'gene'))
    return(GRN_coef.dt)
}