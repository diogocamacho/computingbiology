####################################
# DIFFERENTIAL EXPRESSION via LIMMA
# Diogo M. Camacho
# Wyss Institute 
# Harvard University
# June 2017
#
# Given an expression data and ids for 
# groups to test, runs limma to define 
# differentially expressed genes.
#######################################
limma_dge <- function(expression_data, caseIds=list(), ctrIds)
{
  require(limma)
  if(!is.list(caseIds)) #ie, only one case
  {
    tmp_data <- expression_data[,c(ctrIds,caseIds)]
    
    expdes <- matrix(0,nrow=ncol(tmp_data),2) # experimental design
    expdes[seq(1,length(ctrIds)),1] <- 1
    expdes[seq(length(ctrIds)+1,length(ctrIds)+length(caseIds)),2] <- 1
    colnames(expdes) <- c("control","case")
    
    contmat <- makeContrasts(case-control,levels=expdes) # contrast matrix
    
    limma.fit <- lmFit(tmp_data,design=expdes)
    limma.fit <- contrasts.fit(limma.fit,contmat)
    limma.fit <- eBayes(limma.fit)
    
    res <- topTable(limma.fit,1,number=nrow(tmp_data),sort.by = "none",adjust.method = "fdr")
    return(res)

  } else {

    expdes <- matrix(0,nrow=ncol(expression_data),1)
    expdes[ctrIds,1] <- 1 # controls

    for(i in seq(1,length(caseIds)))
    {
      x <- matrix(0,nrow=ncol(expression_data),1)
      x[caseIds[[i]],1] <- 1
      expdes <- cbind(expdes,x)
    }

    tmp_data <- expression_data[,c(ctrIds,unlist(caseIds))]
    expdes <- expdes[c(ctrIds,unlist(caseIds)),]
    
    colnames(expdes) <- c("control",paste("case",seq(1,length(caseIds)),sep=""))
    myContrasts <- c(paste(colnames(expdes)[2:ncol(expdes)],"-",colnames(expdes)[1],sep=""))
    contmat <- contmat <- eval(as.call(c(as.symbol("makeContrasts"),as.list
                                         (myContrasts),levels=list(expdes)))) # contrast matrix

    limma.fit <- lmFit(tmp_data,design=expdes)
    limma.fit <- contrasts.fit(limma.fit,contmat)
    limma.fit <- eBayes(limma.fit)

    return(limma.fit)
  }
}


limma_voom <- function(count_data,caseIds,ctrIds)
{
  require(limma)
  
  tmp_data <- count_data[,c(ctrIds,caseIds)]
  # samps <- factor(c(rep("control",length(ctrIds)),rep("case",length(caseIds))))
  # expdes <- model.matrix(~0 + samps)
  # colnames(expdes) <- c("control","case")
  
  expdes <- matrix(0,nrow=ncol(tmp_data),2) # experimental design
  expdes[seq(1,length(ctrIds)),1] <- 1
  expdes[seq(length(ctrIds)+1,length(ctrIds)+length(caseIds)),2] <- 1
  colnames(expdes) <- c("control","case")
  
  contmat <- makeContrasts(case-control,levels=expdes) # contrast matrix
  
  V = voom(counts=tmp_data,design=expdes,normalize.method="quantile")
  
  limma.fit <- lmFit(V,design=expdes)
  limma.fit <- contrasts.fit(limma.fit,contmat)
  limma.fit <- eBayes(limma.fit)
  
  res <- topTable(limma.fit,1,number=nrow(tmp_data),sort.by = "none",adjust.method = "fdr")
  
  return(res)
  
}


volcano_dge <- function(limma_res,p_thr,fold_thr)
{
  if(!require(ggplot2)) 
  {
    install.packages("ggplot2")
    library(ggplot2)
  }
  
  if(missing(p_thr)) p_thr <- 0.01
  if(missing(fold_thr)) fold_thr <- 1
  
  ggplot(data=limma_res,aes(x=logFC,y=-log10(adj.P.Val))) + 
    geom_point() + 
    geom_point(data=subset(limma_res,abs(logFC) > fold_thr & adj.P.Val < p_thr),color="red") +
    labs(x="log2(fold change)",y="-log10(fdr p-value)") + 
    xlim(c(-max(abs(limma_res$logFC)),max(abs(limma_res$logFC)))) + 
    theme_bw() + 
    theme(axis.title = element_text(size=20,face="bold"),
          axis.text = element_text(size=12))
  
}

# sorted_dge <- function(limma_res,p_thr,fold_thr)
# {
#   if(!require(ggplot2)) 
#   {
#     install.packages("ggplot2")
#     library(ggplot2)
#   }
#   
#   if(missing(p_thr)) p_thr <- 0.01
#   if(missing(fold_thr)) fold_thr <- 1
#   
#   ggplot() + 
#     geom_point(data=subset(limma_res,adj.P.Val < p_thr),color="red") +
#     geom_point(data=limma_res,aes(y=sort(logFC),x=seq(1,nrow(limma_res)))) + 
#     # labs(x="log2(fold change)",y="-log10(fdr p-value)") + 
#     # xlim(c(-max(abs(limma_res$logFC)),max(abs(limma_res$logFC)))) + 
#     theme_bw() + 
#     theme(axis.title = element_text(size=20,face="bold"),
#           axis.text = element_text(size=12))
#   
# }