.convertToMatrix <- function(mat){
  colNames <- colnames(mat)
  if(inherits(mat, what="DelayedMatrix")){
    # to avoid warning message
    mat <- Matrix::Matrix(mat, ncol=ncol(mat))
  }
  else if(is.matrix(mat)){
    mat <- Matrix::Matrix(mat)
  }
  colnames(mat) <- colNames
  return(mat)
}

.marginMax <- function(mat, margin=c("row", "col")){

  margin <- match.arg(margin, choices=c("col", "row"))
  if(margin=="row"){fun <- MatrixGenerics::rowMaxs}
  else{fun <- MatrixGenerics::colMaxs}

  if(!is(mat, "CsparseMatrix") & !is(mat, "TsparseMatrix")){
    mat <- as.matrix(mat)
  }
  marginMax <- fun(mat)
  return(marginMax)
}

.marginSum <- function(mat, margin=c("row", "col")){

  margin <- match.arg(margin, choices=c("col", "row"))
  if(margin=="row"){fun <- MatrixGenerics::rowSums}
  else{fun <- MatrixGenerics::colSums}

  if(!is(mat, "CsparseMatrix") & !is(mat, "TsparseMatrix")){
    mat <- as.matrix(mat)
  }
  marginMax <- fun(mat)
  return(marginMax)
}

.marginQuant <- function(mat, probs, margin=c("row", "col")){

  margin <- match.arg(margin, choices=c("col", "row"))
  if(margin=="row"){fun <- MatrixGenerics::rowQuantiles}
  else{fun <- MatrixGenerics::colQuantiles}
  if(!is(mat, "CsparseMatrix") & !is(mat, "TsparseMatrix")){
    mat <- as.matrix(mat)
  }
  marginQuant <- t(fun(mat, probs=probs))

  return(marginQuant)
}

.getType <- function(atacFrag) {
  cuts=c(0L,120L,300L,500L)
  atacFrag[,width:=end-start]
  atacFrag[,frag_type:=cut(width, breaks=c(cuts, Inf), labels=TYPENAMES,
                           right=TRUE, include.lowest=TRUE)]
  atacFrag$width <- NULL
  return(atacFrag)
}

#' @author Pierre-Luc Germain
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
                    strandCol="strand", stranded=FALSE, addMetaCols=FALSE){
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

  if(addMetaCols){
    metaCols <- dt[,setdiff(colnames(dt),
                            c(seqCol, startCol, endCol, strandCol,
                              "seqnames", "chr")),with=FALSE]
    mcols(gr) <- metaCols
  }

  return(gr)
}

.getPRCurve <- function(dt,
                        scores="score",
                        labels="cond",
                        models=NULL,
                        posClass="pos",
                        negClass="neg",
                        subSample=FALSE,
                        aggregate=FALSE,
                        seed=42){
  set.seed(seed)

  if(is.null(models)){
    dt$model <- "all"
    models <- "model"
  }

  dt <- copy(dt)
  setnames(dt, scores, "scores")
  setnames(dt, labels, "labels")

  setorder(dt, -scores)
  if(subSample)
  {
    dt <- dt[,.SD[sample(.N, min(1e5,.N))], by=c(models)]
  }

  dt[,tpr:=cumsum(labels==posClass)/sum(labels==posClass), by=c(models)]
  dt[,fpr:=cumsum(labels==negClass)/sum(labels==negClass), by=c(models)]
  dt[,fdr:=cumsum(labels==negClass)/seq_len(.N), by=c(models)]
  dt[,ppv:=cumsum(labels==posClass)/seq_len(.N), by=c(models)]
  dt[,p:=seq_len(.N), by=c(models)]
  dt[,idx:=1:.N, by=c(models)]

  # compute distance to top right
  dt[,dist:=sqrt((1-ppv)^2+(1-tpr)^2)]

  # find sparsification threshold
  dt$is_opt <- FALSE
  dt[which.min(dist), is_opt:= TRUE, by=c(models)]
  dt[, thr:=data.table::first(scores[is_opt == TRUE]), by = c(models)]

  # get precision and recall at optimum
  dt[,tpr_thr:=sum(scores>=thr & labels==posClass)/sum(.SD==posClass),
     by=c(models)]
  dt[,ppv_thr:=sum(scores>=thr & labels==posClass)/sum(scores>=thr),
     by=c(models)]

  dt[,sum_pos:=sum(labels==posClass), by=c(models)]
  dt[,sum_neg:=sum(labels==negClass), by=c(models)]
  dt <- subset(dt, sum_pos>0 & sum_neg>0)

  if(nrow(dt)>0){
    dt[,auc_pr_mod:=PRROC::pr.curve(scores.class0=scores,
                                    weights.class0=as.integer(labels),
                                    curve=FALSE)$auc.integral, # actually needs to be defined whats changing and whats not
       by=c(models)]

    setnames(dt, c("labels"), c(labels))

    if(aggregate){
      dt <- dt[,.(auc_pr_mod=data.table::first(auc_pr_mod),
                  pr=data.table::first(ppv_thr),
                  thr=data.table::first(thr),
                  recall=data.table::first(tpr_thr)), by=c(models)]
    }

    return(dt)
  }
  else
  {
    return(NULL)
  }
}
