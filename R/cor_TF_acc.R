#' Correlate TF expression and Region accessibility
#'
#' This function calculates the correlation between the expression of transcription factors (TFs) and the accessibility of their binding sites across different genomic regions.
#'
#' @param rna.sce A SingleCellExperiment object containing RNA-seq data.
#' @param atac.sce A SingleCellExperiment object containing ATAC-seq data.
#' @param TFs_filt Character vector containing the names of TFs to include in the analysis. Default is NULL.
#' @param motifmatcher.se SingleCellExperiment object containing motif matches.
#' @param motif2gene.dt Data table containing the mapping of motifs to genes.
#' @param correlation_method Method to calculate correlation. Default is "pearson", other options include "spearman" and "kendall".
#' @param remove_motifs character vector of motifs to be manually removed
#' @param cores Number of cores to use for parallel processing.
#'
#' @return A SummarizedExperiment object containing the correlation matrix and the p-value matrix.
#'
#' @import data.table
#' @import dplyr
#' @import SingleCellExperiment
#' @importFrom dplyr bind_cols
#' @importFrom psych corr.test
#' @importFrom tibble column_to_rownames
#' @importFrom parallel detectCores mclapply
#'
#' @export

cor_TF_acc = function(rna.sce, 
                      atac.sce, 
                      TFs_filt = NULL,
                      motifmatcher.se = motifmatcher.se, 
                      motif2gene.dt = motif2gene.dt,
                      correlation_method = "pearson", 
                      remove_motifs = c("T_789"),
                      cores = detectCores()){
    
    ## Filter TFs 
    motifs <- intersect(colnames(motifmatcher.se),motif2gene.dt$motif)
    motifmatcher.se <- motifmatcher.se[,motifs]
    motif2gene.dt <- motif2gene.dt[motif%in%motifs]

    genes <- intersect(toupper(rownames(rna.sce)),motif2gene.dt$gene)
    rna_tf.sce <- rna.sce[genes,]
    rownames(rna_tf.sce) <- toupper(rownames(rna_tf.sce))
    motif2gene.dt <- motif2gene.dt[gene%in%genes]

    # Manually remove some motifs
    if(!is.null(remove_motifs)){
        motif2gene.dt <- motif2gene.dt[!motif %in% remove_motifs]
    }
        
    # Remove duplicated gene-motif pairs
    genes.to.remove <- names(which(table(motif2gene.dt$gene)>1))
    print(sprintf("Removing %d TFs that have duplicated gene-motif pairs:\n%s", length(genes.to.remove), paste(genes.to.remove, collapse=", ")))
    motif2gene.dt <- motif2gene.dt[!gene%in%genes.to.remove]
    rna_tf.sce <- rna_tf.sce[rownames(rna_tf.sce)%in%motif2gene.dt$gene]
    motifmatcher.se <- motifmatcher.se[,colnames(motifmatcher.se)%in%motif2gene.dt$motif]
    stopifnot(table(motif2gene.dt$gene)==1)

    # Sanity checks
    stopifnot(colnames(rna_tf.sce)==colnames(atac.sce))

    TFs <- rownames(rna_tf.sce)

    if(!is.null(TFs_filt)){ 
        TFs = TFs[TFs %in% TFs_filt]
    }
    
    ## Correlate TF-expr & Region-accessibility
    message(sprintf('Correlate TF-expr & Region-accessibility for %s TFs', length(TFs)))

    # a = Sys.time()
    matrix = mclapply(TFs, function(i){
      cat(i, '\n')
      # Get motif name
      motif_i <- motif2gene.dt[gene==i,motif]

      # Get all peaks with motif
      all_peaks_i <- rownames(motifmatcher.se)[which(assay(motifmatcher.se[,motif_i],"motifMatches")==1)]

      # calculate correlations between TF expression and accessibility of motif containing peaks
      corr_output <- psych::corr.test(
        x = t(logcounts(rna_tf.sce[i,])), 
        y = t(assay(atac.sce[all_peaks_i,],"logcounts")), 
        ci = FALSE,
        method = correlation_method
      )

      # create data.table containing cor & pval
      results.dt = data.table(peak = colnames(corr_output$r),
                              cor = round(corr_output$r[1,],3),
                              pval = round(corr_output$p[1,],10)) %>%
        setnames(c('cor', 'pval'), c(paste0('cor.', i), paste0('pval.', i))) %>% 
            .[match(rownames(atac.sce), peak), ] %>% # order peaks by original order
            .[,peak:=NULL]

        return(results.dt)
    }, mc.cores=cores) %>% dplyr::bind_cols(.) %>% # combine all data.tables
        .[,peak := rownames(atac.sce)] %>%  # add peak names
        as.data.frame() %>% 
        tibble::column_to_rownames('peak') %>% as.matrix() # convert to matrix 

    # Extract correlation matrix
    cor.mtx = matrix[,grep('cor', colnames(matrix))]
    colnames(cor.mtx) = gsub('cor.', '', colnames(cor.mtx))

    # Extract Pvalue matrix
    pvalue.mtx = matrix[,grep('pval', colnames(matrix))]
    colnames(pvalue.mtx) = gsub('pval.', '', colnames(cor.mtx))
    
    output = SummarizedExperiment(
          assays = SimpleList("cor" = dropNA(cor.mtx), "pvalue" = dropNA(pvalue.mtx)),
          rowData = rowData(atac.sce)
        )
    return(output)
}