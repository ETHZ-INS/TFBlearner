.getInsertsPos <- function(atacFrag, motifData, stranded, shiftLeft){

  if(stranded){
    specColsMotif <- c("motif_center", "start",
                       "end", "motif_id", "motif_match_id", "strand")}
  else{
    specColsMotif <- c("motif_center", "start",
                       "end", "motif_id", "motif_match_id")}

  if(stranded) specColsFrag <- c("sample", "strand") else specColsFrag <- c("sample")

  # convert to granges for faster overlapping
  motifMarginRanges <- .dtToGr(motifData, startCol="start_margin",
                              endCol="end_margin", seqCol="chr",
                              stranded=stranded)
  atacStartRanges <- .dtToGr(atacFrag, startCol="start", endCol="start",
                            seqCol="chr", stranded=stranded)
  atacEndRanges <- .dtToGr(atacFrag, startCol="end", endCol="end",
                          seqCol="chr", stranded=stranded)

  startHits <- GenomicRanges::findOverlaps(atacStartRanges, motifMarginRanges,
                                           type="within", ignore.strand=TRUE)
  endHits <- GenomicRanges::findOverlaps(atacEndRanges, motifMarginRanges,
                                         type="within", ignore.strand=TRUE)

  # get overlapping insertion sites
  atacStartInserts <- atacFrag[queryHits(startHits),
                               c(specColsFrag, "start"), with=FALSE]
  atacEndInserts <- atacFrag[queryHits(endHits),
                             c(specColsFrag, "end"), with=FALSE]
  setnames(atacStartInserts, "start", "insert")
  if(stranded) setnames(atacStartInserts, "strand", "strand_insert")
  setnames(atacEndInserts, "end", "insert")
  if(stranded) setnames(atacEndInserts, "strand", "strand_insert")

  ai <- cbind(rbindlist(list(atacStartInserts, atacEndInserts)),
              motifData[c(subjectHits(startHits),
                          subjectHits(endHits)), specColsMotif, with=FALSE])

  # count insertions around motif
  ai[,rel_pos:=insert-motif_center]
  ai[,type:=fifelse(insert>=start & insert<=end, 1,0)]
  ai[,ml:=end-start+1, by=motif_id] # motif length

  if(stranded){
    # take strandedness of fragment into account
    if(nrow(ai)>0){
      if(data.table::first(ai$ml) %% 2!=0){
        ai[,rel_pos:=fifelse(strand_insert=="-", -1*rel_pos, rel_pos)]
      }
      else{
        aiMotif <- subset(ai, type==1)
        if(shiftLeft){
          ai[,rel_pos:=fifelse(strand_insert=="-", -1*rel_pos-1, rel_pos)]
          ai[,rel_pos:=fifelse(strand=="-", -1*rel_pos-1, rel_pos)]
        }
        else{
          ai[,rel_pos:=fifelse(strand_insert=="-", -1*rel_pos+1, rel_pos)]
          ai[,rel_pos:=fifelse(strand=="-", -1*rel_pos+1, rel_pos)]
        }
      }}
  }

  return(ai)
}

#' Tn5 insertion counting
#'
#' Counts Tn5 insertions around provided motif-matches.
#' If requested also computes insertion footprint profiles and weighted insertion counts.
#'
#' @name getInsertionProfiles
#' @param atacData Named list of [GenomicRanges::GRanges-class], [data.table::data.table], data.frames or paths to .bed /. bam files
#' containing ATAC-seq fragment coordinates (i.e. chr/seqnames, start, end and optionally a strand column). List names will be used as sample names.
#' If a single object is provided and it contains a column named "sample", insertion counts and profiles will be computed for each sample.
#' @param motifRanges [GenomicRanges::GRanges-class] object containing coordinates of motif-matches.
#' @param margin Margin around motif-matches to consider for computing Tn5 insertion events
#' @param shift If Tn5 insertion bias should be considered (only if strand column is provided).
#' @param calcProfile If insertion footprint profiles should be computed.
#' @param profiles Pre-computed insertion footprint profile to compute weighted insertion counts in case `calcProfile=FALSE`.
#' Needs to contain coordinate (chr/seqnames, start, end) columns and weight column (termed "w").
#' @param symmetric If transcription factor footprint profiles should be symmetric around the motif matches. Only used if `calcProfile=TRUE`.
#' @param stranded If insertion footprint profiles should be computed taking strandedness of fragments into account.
#' @param subSample If fragments should be sub-sampled (to a total of 1e8 fragments per sample) for speed-up.
#' @param BPPARAM Parallel back-end to be used. Passed to [BiocParallel::bpmapply()].
#' @return [data.table::data.table] containing insertion counts within and in margins around motif matches and weighted insertion counts in case
#' an insertion profile is provided or if `calcProfile=TRUE`.
# If `calcProfile=TRUE` also a footprint profile around the motif matches is returned, containing a weight ("w") column corresponding to relative insertion
# frequency at the respective position relative to the motif center ("rel_pos").
#' @import data.table
#' @importFrom GenomicRanges findOverlaps GPos resize GRanges
#' @importClassesFrom GenomicRanges GRanges
#' @export
getInsertionProfiles <- function(atacData,
                                 motifRanges,
                                 margin=200,
                                 shift=FALSE,
                                 calcProfile=TRUE,
                                 profiles=NULL,
                                 symmetric=FALSE,
                                 stranded=FALSE,
                                 subSample=FALSE,
                                 BPPARAM=SerialParam()){

  # prep motif data
  motifData <- .processData(motifRanges, shift=FALSE, readAll=FALSE)
  if(!("motif_id" %in% colnames(motifData))){
    message("Assuming all ranges are of the same type")
    motifData[,motif_id:=1L]
  }

  if(margin>0){
    motifMarginRanges <- as.data.table(GenomicRanges::resize(motifRanges,
                                                             width=2*margin,
                                                             fix="center"))
  }
  else{
    motifMarginRanges <- as.data.table(motifRanges)
  }
  setnames(motifMarginRanges, c("start", "end"), c("start_margin", "end_margin"))
  motifData <- cbind(motifData, motifMarginRanges[,c("start_margin", "end_margin"), with=FALSE])

  # prep ATAC fragment data
  if(is.data.table(atacData)) atacData <- list(atacData)
  atacFrag <- lapply(atacData, .processData, shift=shift, subSample=subSample)
  if(!("sample" %in% colnames(atacFrag[[1]]))){
    names(atacFrag) <- names(atacData)
    atacFrag <- rbindlist(atacFrag, idcol="sample") 
  }
  else{
    atacFrag <- rbindlist(atacFrag)
  }

  commonChr <- intersect(unique(motifData$chr),unique(atacFrag$chr))
  chrLevels <- commonChr

  atacFrag <- subset(atacFrag, chr %in% chrLevels)
  motifData <- subset(motifData, chr %in% chrLevels)

  motifLevels <- unique(motifData$motif_id)

  # convert to factors (memory usage)
  motifData[,chr:=as.integer(factor(chr, levels=chrLevels, ordered=TRUE))]
  motifData[,motif_id:=factor(motif_id, levels=motifLevels, ordered=TRUE)]

  # determine motif center
  motifData[,motif_center:=floor((end-start)/2)+start]

  motifData[,end_margin:=fifelse(end_margin-motif_center<margin,end_margin+1,end_margin)]
  motifData[,start_margin:=fifelse(start_margin-motif_center>-margin,start_margin-1,start_margin)]

  distEnd <- motifData$end[1] - motifData$motif_center[1]
  distStart <- motifData$motif_center[1] - motifData$start[1]
  if(distStart>distEnd){
    shiftLeft <- TRUE}
  else{
    shiftLeft <- FALSE
  }

  atacFrag[,chr:=as.integer(factor(chr, levels=chrLevels, ordered=TRUE))]
  medZero <- function(x, len){ median(c(rep(0,max(0,len-length(x))),x)) }
  nSamples <- length(unique(atacFrag$sample))

  setorder(motifData, chr)
  setorder(atacFrag, chr)
  motifData[,motif_match_id:=1:nrow(motifData)]

  motifData <- split(motifData, by="chr")
  atacFrag <- split(atacFrag, by="chr")

  if(calcProfile){
   atacProfiles <- BiocParallel::bpmapply(function(md, af, stranded, shiftLeft){
      atacInserts <- .getInsertsPos(af, md, stranded, shiftLeft)
      atacProfile <- atacInserts[,.(pos_count_global=.N),
                                 by=.(ml, rel_pos, sample, motif_id, type)]
      return(atacProfile)
   }, motifData, atacFrag, MoreArgs=list(stranded=stranded, shiftLeft=shiftLeft), 
    SIMPLIFY=FALSE, BPPARAM=BPPARAM)
    atacProfiles <- rbindlist(atacProfiles, idcol="seqnames")

    atacProfiles <- atacProfiles[,.(pos_count_global=sum(pos_count_global)),
                                 by=.(rel_pos, motif_id, type, ml)]
    # get inserts within motif
    atacProfilesMotif <- subset(atacProfiles, type==1)
    atacProfilesMotif[,med_pos_count_global:=medZero(pos_count_global,
                                                     data.table::first(ml)),
                      by=.(motif_id)]
    atacProfilesMotif <- atacProfilesMotif[,
   .(pos_count_global=(pos_count_global+data.table::first(med_pos_count_global))/2),
                                           by=.(rel_pos, motif_id, type)]

    atacProfilesMargin <- subset(atacProfiles, type==0)

    # fill non covered positions
    atacProfiles <- rbind(atacProfilesMargin, atacProfilesMotif, fill=TRUE)
    allPos <- data.table(expand.grid(motifLevels, seq(-margin,margin)))

    colnames(allPos) <- c("motif_id", "rel_pos")
    allPos$motif_id <- factor(allPos$motif_id, levels=motifLevels, ordered=TRUE)
    allPos$pos_count_global <- 0

    colsProfile <- c("rel_pos", "motif_id", "pos_count_global")
    atacProfiles[,rel_pos:=as.integer(rel_pos)]

    atacProfiles <- rbind(atacProfiles[, colsProfile, with=FALSE],
                          allPos[!atacProfiles, on=c("rel_pos", "motif_id")])

    # calculate weights
    setorder(atacProfiles, motif_id, rel_pos)
    atacProfiles[,w:=smooth(pos_count_global, twiceit=TRUE), by=motif_id]
    if(symmetric) atacProfiles[,w:=rev(w)+w, by=motif_id]
    atacProfiles[,w:=length(w)*w/sum(w), by=motif_id]
    atacProfiles[,w_count:=pos_count_global/max(pos_count_global), by=motif_id]
    atacProfiles[,w_count_smooth:=smooth(w_count), by=motif_id]
  }
  else{
    atacProfiles <- profiles
  }
  
  # get match scores
  motifScores <- BiocParallel::bpmapply(function(md,af,
                                                 stranded,
                                                 profiles,
                                                 shiftLeft){

    atacInserts <- .getInsertsPos(af, md, stranded, shiftLeft)

    if(!is.null(profiles)){
      atacInserts <- atacInserts[,.(pos_count=.N),
                                 by=.(motif_match_id, motif_id, sample, 
                                      rel_pos, type)]
      atacInserts <- merge(atacInserts,
                           profiles[, c("rel_pos", "motif_id", "w"),with=FALSE],
                           by.x=c("motif_id","rel_pos"),
                           by.y=c("motif_id","rel_pos"), all.x=TRUE, all.y=FALSE)
      atacInserts[,score:=w*pos_count]
      atacInsertSum <- atacInserts[,.(score=sum(score),
                                      tot_count=sum(pos_count)),
                                   by=.(motif_match_id, motif_id, sample, type)]
    }
    else{
      atacInsertSum <- atacInserts[,.(tot_count=.N),
                                   by=.(motif_match_id, motif_id, sample, type)]
    }
    return(atacInsertSum)
  }, motifData, atacFrag,
  MoreArgs=list(stranded=stranded,
                profiles=atacProfiles,
                shiftLeft=shiftLeft),
  SIMPLIFY=FALSE, 
  BPPARAM=BPPARAM)

  motifScores <- rbindlist(motifScores)

  motifData <- rbindlist(motifData)
  motifData[,chr:=chrLevels[chr]]
  motifScores <- cbind(motifScores,
                       motifData[motifScores$motif_match_id,
                                 c("start", "end", "chr"), with=FALSE])

  if("tot_count" %in% colnames(motifScores)){
  setnames(motifScores, c("tot_count"), c("insert_counts"))}
  if("score" %in% colnames(motifScores)){
    setnames(motifScores, c("score"), c("weighted_insert_counts"))}

  return(list(motifScores=motifScores, profile=atacProfiles))
}
