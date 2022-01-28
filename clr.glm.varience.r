clr.glm.varience.local <- function(tdata, rdata, clr_output, top_tfs, genesName, tfsName, intercept) {
  genes <- scan(genesName, character())
  tfs <- scan(tfsName, character())
  #tdata dimensions: genes by samples
  #rdata dimensions: tfs by samples
  #clr_output dimensions: tfs by genes
  cat(as.character(Sys.time()),"\n")
  #b.adj dimensions: tfs by genes
  B.adj <- matrix(0,nrow=length(tfs),ncol=length(genes))
  #filter out tf by tf matches
  for (gene_index in 1:length(genes)) {
    for (tf_index in 1:length(tfs)) {
      if (genes[gene_index] == tfs[tf_index]) {
        clr_output[tf_index, gene_index] = 0;
      }
    }
  }
  cat("Working on Gene:")
  for(gene_index in 1:length(genes)) {#for each gene
    cat(gene_index,", ")
    if (gene_index == 84) {
      next
    }
    #run glm
    glm_output <- bestglm(data.frame(cbind(t(rdata[order(clr_output[,gene_index], decreasing = T)[1:top_tfs],]), tdata[gene_index,])), intercept = intercept)[["BestModel"]]$coefficients
    #put output into output matrix
    if (intercept == T) {
      start <- 2
    } else {
      start <- 1
    }
    if (length(glm_output) > 1) {
      for (coeff_index in start:length(glm_output)) {
        B.adj[order(clr_output[, gene_index], decreasing = T)[as.numeric(substring(names(glm_output)[coeff_index], 2))], gene_index] <- glm_output[coeff_index]
      }
    }
  }
  #transpose output matrix and multiply by tf levels to get genes by samples then subtract from actual
  #errorScores dimensions: tfs by genes
  errorScores <- matrix(0,nrow=dim(rdata)[1],ncol=dim(tdata)[1])#create vector
  errorVector <- 1:length(genes)
  cat("Working on Gene:")
  for (gene_index in 1:dim(tdata)[1]) {
    cat(gene_index,", ")
    totalError <- sum((tdata [gene_index,] - t(B.adj [,gene_index]) %*% rdata)^2)
    errorVector [gene_index] <- totalError
    for(tf_index in 1:dim(rdata)[1]) {
      tfRemovedCoefficients <- t(B.adj [,gene_index])
      tfRemovedCoefficients [tf_index] <- 0#create matrix without 1 tf
      errorScores[tf_index, gene_index] <- 1 - totalError/sum((tdata [gene_index,] - tfRemovedCoefficients %*% rdata)^2)#calculate score
    }
  }
  cat("\n")
  # 	##return results
  errorScores
}
