#' Predict TF binding sites and update motifmatchr results using virtual ChIP-seq library
#'
#' @param atac.sce A SingleCellExperiment object containing ATAC-seq data.
#' @param tf2peak_cor.se SummarizedExperiment object containing correlation scores between TFs and peaks (output of cor_TF_acc())
#' @param motifmatcher.se SingleCellExperiment object containing motif matches.
#' @param motif2gene.dt Data table containing the mapping of motifs to genes.
#' @param min_number_peaks minimum number of peaks required to create a virtual ChIP-seq library for a TF
#' @param TFs_filt Character vector containing the names of TFs to include in the analysis. Default is NULL.
#' @param remove_motifs character vector of motifs to be manually removed
#' @param cores integer indicating number of CPU cores to be used for parallel processing
#'
#' @return a list containing a data.table with virtual_chip.dt and a SingleCellExperiment object containing updated motif match scores
#'
#' @import data.table
#' @import dplyr
#' @import SingleCellExperiment
#' @importFrom Matrix Matrix
#' @importFrom parallel detectCores mclapply
#'
#' @export

silico_chip = function(atac.sce = atac.sce, 
                       tf2peak_cor.se = tf2peak_cor.se,
                       motifmatcher.se = motifmatcher.se,
                       motif2gene.dt = motif2gene.dt, 
                       min_number_peaks = 50,
                       TFs_filt = NULL, 
                       remove_motifs = NULL,
                       cores = detectCores()){
    
    ## Subset peaks 
    peaks <- intersect(rownames(motifmatcher.se),rownames(tf2peak_cor.se))
    print(sprintf("Number of peaks: %s",length(peaks)))

    tf2peak_cor.se <- tf2peak_cor.se[peaks,]
    atac.sce <- atac.sce[peaks,]
    motifmatcher.se <- motifmatcher.se[peaks,]
    
    # Manually remove some motifs
    if(!is.null(remove_motifs)){
        motif2gene.dt <- motif2gene.dt[!motif %in% remove_motifs]
    }

    ## Rename TFs 
    motifs <- intersect(colnames(motifmatcher.se),motif2gene.dt$motif)
    TFs <- intersect(colnames(tf2peak_cor.se),motif2gene.dt$gene)
    if(!is.null(TFs_filt)){ TFs = TFs[TFs %in% TFs_filt]}

    motif2gene_filt.dt <- motif2gene.dt[motif%in%motifs & gene%in%TFs]
    motifs <- motif2gene_filt.dt$motif
    TFs <- motif2gene_filt.dt$gene

    tmp <- TFs; names(tmp) <- motifs

    stopifnot(motif2gene_filt.dt$motif%in%colnames(motifmatcher.se))
    stopifnot(motif2gene_filt.dt$gene%in%colnames(tf2peak_cor.se))

    motifmatcher.se <- motifmatcher.se[,motifs]
    colnames(motifmatcher.se) <- tmp[colnames(motifmatcher.se)]
    tf2peak_cor.se <- tf2peak_cor.se[,TFs]
    stopifnot(colnames(motifmatcher.se)==colnames(tf2peak_cor.se))

    ## Prepare data 
    tf2peak_cor.mtx <- assay(tf2peak_cor.se,"cor")
    motifmatcher.mtx <- assay(motifmatcher.se,"motifScores")
    atac.mtx <- assay(atac.sce[peaks,],"logcounts") %>% round(3)
    
    ######################################
    ## Create virtual chip-seq library ##
    ######################################

    print("Predicting TF binding sites...")
    stopifnot(!duplicated(TFs))
    print(sprintf("Number of TFs: %s",length(TFs)))

    virtual_chip.dt <- mclapply(TFs, function(i) {
      print(i)

      peaks <- names(which(abs(tf2peak_cor.mtx[,i])>0)) # we only consider chromatin activators # The 'abs' actually makes it so that all non-zero are kept (== all)

      if (length(peaks)>=min_number_peaks) {

        # calculate accessibility score
        max_accessibility_score <- apply(atac.mtx[peaks,],1,max) %>% round(2)

        # calculate correlation score
        correlation_score <- tf2peak_cor.mtx[peaks,i] %>% round(2)
        correlation_score[correlation_score==0] <- NA

        # calculate motif score
        motif_score <- motifmatcher.mtx[peaks,i]
        motif_score <- round(motif_score/max(motif_score),2)

        # calculate motif counts
        predicted_score <- correlation_score * minmax.normalisation(max_accessibility_score * motif_score)

        tmp <- data.table(
          peak = peaks, 
          correlation_score = correlation_score,
          max_accessibility_score = max_accessibility_score,
          motif_score = motif_score,
          score = round(predicted_score,2)
        ) %>% sort.abs("score") %>% 
          .[,c("peak","score","correlation_score","max_accessibility_score","motif_score")]

        to_return.dt <- tmp[!is.na(score),c("peak","score")]  %>% .[,tf:=i]
        return(to_return.dt)
      }
    }, mc.cores=cores) %>% data.table::rbindlist(.)
    
    output = list()
    output$virtual_chip.dt = virtual_chip.dt
    
    # Create Virtual ChIP-seq matrix
    virtual_chip.mtx <- virtual_chip.dt %>% 
          .[,peak:=factor(peak,levels=rownames(motifmatcher.se))] %>%
          data.table::dcast(peak~tf, value.var="score", fill=0, drop=F) %>%
          matrix.please %>% Matrix::Matrix(.)
    
    ## Update motifmatchr results using the virtual ChIP-seq library
    print("Updating motifmatchr results using the virtual ChIP-seq library...")

    motifmatcher_chip.se <- motifmatcher.se[,colnames(virtual_chip.mtx)]
    # reset motif matches
    stopifnot(rownames(virtual_chip.mtx)==rownames(motifmatcher_chip.se))
    stopifnot(colnames(virtual_chip.mtx)==colnames(motifmatcher_chip.se))
    assay(motifmatcher_chip.se,"VirtualChipScores") <- virtual_chip.mtx
    
    output$motifmatcher_chip.se = motifmatcher_chip.se
    
    return(output)
}