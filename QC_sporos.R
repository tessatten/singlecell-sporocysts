#date: 10-March-2020
#Project: Basic QC for whole cells

################
#load libraries#
################
library(Seurat)
library(reticulate)
#conda_install("r-reticulate", "scrublet")

#######
#Usage#
#######

#filtered = wholeQC("your project directory")

#it will save the original seurat ("QC_object_pass_fail_all_cells.RDS") with PASS or Fail based on the QC metrics below
#it will also save the filtered seurat ("QC_object_pass_only_cells.RDS") with PASS or Fail based on the QC metrics below
#it will return the filtered seurat 




################
#Basic pipeline#
################

#' Basic QC and filtering functions.  You will need reticulate and scrublet installed. 
#' Initial filtering and metadata
#' Calculates metadata, marks cells based on QC metrics.  Assumes Seurat V3 and assumes scrublet
#' @param projectDir directories to load. Must contain the matrix, barcode and features as given by Cell Ranger.
#' @param sGenes Genes used to identify S phase genes.  If NULL, will be loaded from Seurat.
#' @param g2mGenes Genes used to identify G2M phase genes.  If NULL, will be loaded from Seurat.
#' @param maxMT Maximum percentage expression from mitochondrial genes allow per cell.
#' @param minGenes Minimum number of genes for a cell to express and pass QC.
#' @param minUMIs Minimum number of UMIs for a cell to have and pass QC.
#' @param scrubScoreMax Mark a cell as a doublet if its scrublet score exceeds this value
#' @return A filtered Seurat object after QC filters.
wholeQC = function(projectDir,sGenes=NULL,g2mGenes=NULL,maxMT=25,minGenes=200,minUMIs=300,scrubScoreMax = 0.2,numPCs=50,doPlot=TRUE,...){
  if(is.null(projectDir)){
    message("Data directories not named, please name")
  }
  srat = Read10X(projectDir,...)
  srat = CreateSeuratObject(srat)

  # if(is.null(sGenes) | is.null(g2mGenes)){
  #   data('cc.genes.updated.2019',package='Seurat')
  # }
  # if(is.null(sGenes))
  #   sGenes = cc.genes.updated.2019$s.genes
  # if(is.null(g2mGenes))
  #   g2mGenes = cc.genes.updated.2019$g2m.genes
  #Get extra meta-data things
  srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "Smp-9")
  mito.genes <- grep("^MT-", rownames(srat@assays$RNA@data), value = T)

  #Get cell cycle
  srat = NormalizeData(srat,verbose=FALSE)
  #The seeed parameter magic is needed because the Seurat authors are dicks
  # srat = CellCycleScoring(srat, s.features=sGenes,g2m.features=g2mGenes,seed=sample(1e9,1))
  #Call doublets
  srat@meta.data$scrubScore=NA
  srat@meta.data$scrubCall=NA
  srat@meta.data$scrubSim1=NA
  srat@meta.data$scrubSim2=NA
  for(i in seq_along(projectDir)){
    scrubScores = runScrublet(projectDir[i],nPCs=numPCs)
    #Save the relevant information in meta.data
    m = match(scrubScores$barcode,rownames(srat@meta.data))
    s = !is.na(m)
    srat@meta.data$scrubScore[m[s]] = scrubScores$score[s]
    srat@meta.data$scrubCall[m[s]] = scrubScores$call[s]
    srat@meta.data$scrubSim1[m[s]] = scrubScores$simScore1[s]
    srat@meta.data$scrubSim2[m[s]] = scrubScores$simScore2[s]
  }
  #Do basic processing and clustering. This is so that we can plot afterwards but the we will strip the clustering in the filtered object
  srat = FindVariableFeatures(srat,verbose=FALSE)
  srat = ScaleData(srat,verbose=FALSE)
  srat = RunPCA(srat,npcs=numPCs,approx=FALSE,verbose=FALSE)
  srat = FindNeighbors(srat,dims=seq(numPCs),verbose=FALSE)
  srat = FindClusters(srat,res=1,verbose=FALSE)
  srat = RunUMAP(srat,dims=seq(numPCs),verbose=FALSE)

  #Set the basic filters.  Default to TRUE if NA
  srat@meta.data$PASS_MT = !(srat@meta.data$percent.mt >= maxMT)
  srat@meta.data$PASS_nGenes = !(srat@meta.data$nFeature_RNA <= minGenes)
  srat@meta.data$PASS_nCounts = !(srat@meta.data$nCount_RNA <= minUMIs)
  srat@meta.data$PASS_doublet = !(!is.na(srat@meta.data$scrubCall) & (srat@meta.data$scrubCall | srat@meta.data$scrubScore > scrubScoreMax))

  #Work out the final filter
  srat@meta.data$PASS = with(srat@meta.data,PASS_MT & PASS_nGenes & PASS_nCounts & PASS_doublet)

  #Make a reason for fail variable
  nBasicFail = with(srat@meta.data,4 - (PASS_MT + PASS_nGenes + PASS_nCounts + PASS_doublet))
  tmp = rep('Pass',length(nBasicFail))
  tmp[!srat@meta.data$PASS_MT] = '#highMT'
  tmp[!srat@meta.data$PASS_nGenes] = '#Lowgenes'
  tmp[!srat@meta.data$PASS_nCounts] = '#Lowcounts'
  tmp[!srat@meta.data$PASS_doublet] = '#Scrubletdoublet'
  tmp[nBasicFail>1] = 'multiple'

  srat@meta.data$reasonForFail = tmp


  #################
  #Make some Plots#
  #################

    #Plot distributions of basic things
    df = srat@meta.data
    pdf(file.path(projectDir,'QC_features.pdf'))
    plot(density(log10(df$nFeature_RNA)),main='#Genes',xlab='log10(No Genes)')
    abline(v=log10(minGenes),col='red')
    print(plot)
    dev.off()

    pdf(file.path(projectDir,'QC_umis.pdf'))
    plot(density(log10(df$nCount_RNA)),main='#UMIs',xlab='log10(No UMIs)')
    abline(v=log10(minUMIs),col='red')
    print(plot)
    dev.off()

    pdf(file.path(projectDir,'QC_scrublet.pdf'))
    plot(density(df$scrubScore),main='Scrublet Score',xlab='Scrub Score')
    print(plot)
    dev.off()
    
    #Standard QC things
    pdf(file.path(projectDir,'QC_mt_nFeature_nCount.pdf'))
    plot(FeaturePlot(srat,c('percent.mt','nFeature_RNA','nCount_RNA')))
    print(plot)
    dev.off()

    ########
    #Filter#
    ########
  w = which(srat@meta.data$PASS)

  #Make the new version, store the old one in it
  sratOld = srat
  # save the object without filtering
  saveRDS(sratOld, file.path(projectDir,'QC_object_pass_fail_all_cells.RDS'))
  #Restart without them
  mDat = srat@meta.data[w,]
  #Drop old clustering
  mDat = mDat[,!colnames(mDat) %in% c('seurat_clusters',grep('^RNA_snn_res',colnames(mDat),value=TRUE))]
  srat = CreateSeuratObject(srat@assays$RNA@counts[,w])
  #Merge in metadata
  srat@meta.data = cbind(srat@meta.data,mDat[,!(colnames(mDat) %in% colnames(srat@meta.data)),drop=FALSE])
  srat = NormalizeData(srat,verbose=FALSE)
  srat = FindVariableFeatures(srat,verbose=FALSE)
  #stats on the number of cells and genes prior QC
  srat@misc$preQC = sratOld
  # save the filtered object after filtering
  saveRDS(srat, file.path(projectDir,'QC_object_pass_only_cells.RDS'))
  return(srat)
}

#' Run scrublet
#Installing scrublet
# conda_install("r-reticulate", "scrublet")

#' @param dat10X Directory containing 10X matrix and genes files.
#' @param nPCs Number of PCs to use.
#' @return A data.frame with the scrublet results
runScrublet = function(dat10X,nPCs){
  fNom = tempfile()
  cmd = sprintf(
  "import scrublet as scr
import scipy.io
import numpy as np
import os
import gzip
os.chdir('%s')
input_dir = os.path.expanduser('%s')
output_file = '%s'
try:
  counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx').T.tocsc()
  genes = np.array(scr.load_genes(input_dir + '/genes.tsv', delimiter='\\t', column=1))
except FileNotFoundError:
  counts_matrix = scipy.io.mmread(input_dir + '/matrix.mtx.gz').T.tocsc()
  #Copy relevant part of scr.load_genes
  gene_list = []
  gene_dict = {}
  with gzip.open(input_dir + '/features.tsv.gz','rt') as f:
    for l in f:
      gene = l.strip('\\n').split('\\t')[1]
      if gene in gene_dict:
        gene_dict[gene] += 1
        gene_list.append(gene + '__' + str(gene_dict[gene]))
        if gene_dict[gene] == 2:
          i = gene_list.index(gene)
          gene_list[i] = gene + '__1'
      else:
        gene_dict[gene] = 1
        gene_list.append(gene)
  genes = np.array(gene_list)
scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2,
                                                        min_cells=3,
                                                        min_gene_variability_pctl=85,
                                                        n_prin_comps=%d)
with open(output_file,'w') as f:
  #Header row
  f.write('score\\tcall\\tsimScore1\\tsimScore2\\n')
  #As nSim = 2*nObs, we can split into two columns
  offset = len(doublet_scores)
  f.write('\\n'.join([str(x)+'\\t'+str(predicted_doublets[i]).upper()+'\\t'+str(scrub.doublet_scores_sim_[i])+'\\t'+str(scrub.doublet_scores_sim_[offset+i]) for i,x in enumerate(doublet_scores)]))
  f.close()",getwd(),dat10X,fNom,nPCs)
  system(paste0('python -c "',cmd,'"'))
  out = read.table(fNom,header=TRUE,sep='\t')
  if(file.exists(file.path(dat10X,'barcodes.tsv'))){
    tmp = read.table(file.path(dat10X,'barcodes.tsv'),header=FALSE)
  }else{
    tmp = read.table(file.path(dat10X,'barcodes.tsv.gz'),header=FALSE)
  }
  out$barcode = tmp[,1]
  return(out)
}
