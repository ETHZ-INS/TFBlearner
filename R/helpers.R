.getType <- function(atacFrag, cuts=c("nucleosome_free"=0,
                                      "mononucleosome"=120,
                                      "dinucleosome"=300,
                                      "multinucleosome"=500),
                     label=FALSE) {
  atacFrag[,width:=end-start]
  if(label){
    labels <- names(cuts)
  }
  else{
    labels <- FALSE
  }
  atacFrag[,frag_type:=cut(width, breaks=c(cuts, Inf), labels=labels,
                           right=TRUE, include.lowest=TRUE)]
  return(atacFrag)
}

.getNonRedundantMotifs <- function(tfNames,
                                   format=c("PFMatrix","universal","PWMatrix"),
                                   species=c("Hsapiens","Mmusculus")){
  species <- match.arg(species)
  motifs <- MotifDb::query(MotifDb::MotifDb,
                           andStrings=c(species,"HOCOMOCO"),
                           orStrings=tfNames)
  pat <- paste0("^",species,"-HOCOMOCOv1[0-1]-|_HUMAN.+|_MOUSE.+|core-[A-D]-|secondary-[A-D]-")
  modf <- data.frame(row.names=names(motifs),
                     TF=gsub(pat,"",names(motifs)),
                     grade=gsub(".+\\.","",names(motifs)))
  modf <- modf[order(modf$TF,-as.numeric(grepl("HOCOMOCOv11",row.names(modf))),modf$grade),]
  modf <- modf[!duplicated(modf$TF),]
  motifs <- motifs[row.names(modf)]
  switch(match.arg(format),
         universal=setNames(universalmotif::convert_motifs(motifs), modf$TF),
         PFMatrix=do.call(TFBSTools::PFMatrixList, setNames(
           universalmotif::convert_motifs(motifs, class="TFBSTools-PFMatrix"),
           modf$TF)),
         PWMatrix=do.call(TFBSTools::PWMatrixList,
                          setNames(universalmotif::convert_motifs(motifs,
                                                                  class="TFBSTools-PWMatrix"),
                                   modf$TF))
  )
}

.dtToGr <- function(dt, seqCol="seqnames", startCol="start", endCol="end",
                    strandCol="strand", stranded=FALSE){
  dt <- copy(dt)
  setnames(dt, seqCol, "seqnames", skip_absent = TRUE)

  if(stranded) strand <- dt[[strandCol]] else strand <- NULL

  gr <- GRanges(seqnames=dt[["seqnames"]],
                strand=strand,
                ranges=IRanges(start=dt[[startCol]], end=dt[[endCol]]))
  if(startCol==endCol)
  {
    gr <- GPos(seqnames=dt[["seqnames"]],
               strand=strand,
               pos=dt[[startCol]])
  }
  return(gr)
}
