#' Run ChromVAR ChIP-seq
#'
#' This function performs ChromVAR on the subset of in silico ChIP-seq peaks. It takes in an atac.sce object and a motifmatcher_chip.se object and outputs the chromVAR scores for in silico ChIP-seq peaks each motif.
#'
#' @param atac.sce A SingleCellExperiment object containing the ATAC-seq data
#' @param motifmatcher_chip.se A MotifMatchR object (output of silico_chip function())
#' @param assay Character string specifying which assay of motifmatcher_chip.se to use. Default is 'VirtualChipScores'
#' @param background A SummarizedExperiment object containing the background peaks. If NULL, the background peaks are calculated using the ATAC-seq data and genome annotation provided in the genome argument. 
#' @param genome BSgenome object used to calculate background peaks when background = NULL or method = 'chromVAR'.
#' @param positive_only Logical indicating whether to filter out negative TF binding values. Default is TRUE.
#' @param min_chip_score Minimum ChIP-seq score to use for filtering motifs. Default is 0.15.
#' @param min_number_peaks Minimum number of peaks to include for each TF. Default is 50.
#' @param TFs_filt Character vector containing the names of TFs to include in the analysis. Default is NULL.
#' @param test Logical indicating whether to use a test set of TFs for the analysis. Default is FALSE.
#' @param method Character string indicating whether to use the ArchR or ChromVAR implementation for the ChIP-seq analysis. Default is 'ArchR'.
#' @param cores Integer specifying the number of cores to use for parallel processing. Default is the number of cores available.
#' @import data.table
#' @import dplyr
#' @import SingleCellExperiment
#' @importFrom GenomicRanges GRanges
#' @importFrom chromVAR addGCBias getBackgroundPeaks
#' @importFrom BiocParallel register MulticoreParam
#' 
#' @return A SummarizedExperiment object with a VirtualChipScores assay containing the ChromVAR scores for each motif.
#' @export

chromVAR_chip = function(atac.sce = atac.sce, 
                           motifmatcher_chip.se = motifmatcher_chip.se,
                           assay = 'VirtualChipScores',
                           background = bgdPeaks.se,
                           genome = BSgenome.Mmusculus.UCSC.mm10, # Only needed if background = NULL or method = 'chromVAR'
                           positive_only = TRUE, 
                           min_chip_score = 0.15,
                           min_number_peaks = 50,
                           TFs_filt = NULL,
                           test = FALSE,
                           method = 'ArchR', # Option 'ArchR' or 'ChromVAR'
                           cores = detectCores()){
    
    if(method == 'ChromVAR'){
        background = NULL
        cat('Background peaks recalculated for method = ChromVAR \n')
        if(is.null(genome)){
            cat('Please provide genome \n')
            stop()
        }
    }
    
    if(method == 'ArchR'){
        archr_installed = requireNamespace("ArchR", quietly = TRUE)
        if(!archr_installed){ 
            cat('ArchR not installed, install ArchR or use method = "ChromVAR"')
            stop()
    }
    
    if (positive_only){
      print(sprintf("Number of matches before filtering negative TF binding values: %d",sum(assay(motifmatcher_chip.se,"motifMatches"))))
      assay(motifmatcher_chip.se,"motifMatches")[assay(motifmatcher_chip.se, assay)<0] <- F
      print(sprintf("Number of matches after filtering negative TF binding values: %d",sum(assay(motifmatcher_chip.se,"motifMatches"))))
    }

    print(sprintf("Number of matches before filtering based on minimum ChIP-seq score: %d",sum(assay(motifmatcher_chip.se,"motifMatches"))))
    assay(motifmatcher_chip.se,"motifMatches")[abs(assay(motifmatcher_chip.se, assay))<=min_chip_score] <- F
    print(sprintf("Number of matches after filtering based on minimum ChIP-seq score: %d",sum(assay(motifmatcher_chip.se,"motifMatches"))))

    assays(motifmatcher_chip.se) <- assays(motifmatcher_chip.se)["motifMatches"]
    
    
    ## Filter TFs 
    # Filter TFs with too few peaks
    TFs <- which(colSums(assay(motifmatcher_chip.se,"motifMatches")) >= min_number_peaks) %>% names
    TFs.removed <- which(colSums(assay(motifmatcher_chip.se,"motifMatches")) < min_number_peaks) %>% names
    
    # Subset to TFs of interest
    if(!is.null(TFs_filt)){ 
        TFs = TFs[TFs %in% TFs_filt]
    }
    if(length(TFs.removed>0)){
        cat(sprintf("%s TFs removed because they don't have enough binding sites: %s \n", length(TFs.removed), paste(TFs.removed, collapse=" ")))
    }
    
    if(test){
      TFs <- c("FOXA2","MIXL1","GATA1","EOMES","BCL11B","DLX2","FOXC1")
    }
    
    motifmatcher_chip.se <- motifmatcher_chip.se[,TFs]    
    stopifnot(rownames(atac.sce)==rownames(motifmatcher.se))  
    
    if(!is.null(background)){
        ## Load background peaks
        bgdPeaks.se = background
        tmp <- rowRanges(bgdPeaks.se)
        rownames(bgdPeaks.se) <- sprintf("%s:%s-%s",seqnames(tmp), start(tmp), end(tmp))
        bgdPeaks.se <- bgdPeaks.se[rownames(atac.sce),] # subset to same peaks as in filtered atac.sce
        bg = assay(bgdPeaks.se)
    } else{
        ## Determine background peaks
        if(is.null(genome)){
            cat('Please provide genome \n')
            stop()
        }
    
        ## Filter non-accessible peaks
        peaks = rowSums(assay(atac.sce, 'counts'))>0
        atac.sce = atac.sce[peaks,]
        motifmatcher_chip.se = motifmatcher_chip.se[peaks,]   
        
        ## Calculate background peaks
        cat('Calculating background peaks \n')
        peaks = rownames(atac.sce_filt)
        gr <- GRanges(
            seqnames = Rle(peaks %>% str_split(':') %>% map_chr(1)),
            ranges = IRanges(start = as.numeric(peaks %>% str_split(':') %>% map_chr(2) %>% str_split('-') %>% map_chr(1)), 
                             end = as.numeric(peaks %>% str_split(':') %>% map_chr(2) %>% str_split('-') %>% map_chr(2))),
            strand = Rle(rep('*', length(peaks))))
        atac.rse = as(atac.sce_filt, 'RangedSummarizedExperiment')
        atac.rse@rowRanges = gr
        
        ## Adding background peaks
        atac.rse <- chromVAR::addGCBias(atac.rse, genome = genome)
        bg <- chromVAR::getBackgroundPeaks(object = atac.rse)
    }

    
    ## ChromVAR ChIP-seq
    if(method == 'ArchR'){
        archr_installed = requireNamespace("ggrastr", quietly = TRUE)
        if(!archr_installed){ 
            cat('ArchR not installed, install ArchR or use method = "ChromVAR"')
            stop()
        }
        cat('Running ChromVAR-ChIP-seq (ArchR implementation) \n')
        
        featureDF <- data.frame(
          rowSums = rowSums(assay(atac.sce))#, 
       #   start = rowData(atac.sce)$start, # -> Check if these can be left out
       #   end = rowData(atac.sce)$end
        )

        # Compute deviations
        chromvar_deviations.se <- ArchR:::.customDeviations(
          countsMatrix = assay(atac.sce),
          annotationsMatrix = as(assay(motifmatcher_chip.se),"dgCMatrix"),
          backgroudPeaks = bg,
          expectation = featureDF$rowSums/sum(featureDF$rowSums),
          prefix = "",
          out = c("deviations", "z"),
          threads = cores,
          verbose = TRUE
        )
    } else if(method == 'ChromVAR'){
        cat('Running ChromVAR-ChIP-seq (ChromVAR implementation) \n')
        #stop()
        register(MulticoreParam(cores))
        assayNames(atac.sce) <- "counts"

        # Compute deviations
        chromvar_deviations_chromvar.se <- chromVAR::computeDeviations(
          object = atac.sce,
          annotations = assay(motifmatcher_chip.se),
          background_peaks = bg
        )
        
    } else{
        print('Choose either option "ArchR" or "ChromVAR"')
    }
        
    return(chromvar_deviations.se)
}
