#!/staging/biology/ls807terra/0_Programs/anaconda3/envs/RNAseq_quantTERRA/bin/Rscript

# ----------------------------------------------------------------------------- #
# Options
pacman::p_load("optparse")

option_list = list(
  make_option(c("-i", "--counts"), type="character", default=NULL,
              help="Please Enter a count table.", metavar="COUNTS"),
  make_option(c("-a", "--annotation"), type="character", default=NULL,
              help="Please specify a annotation table in csv or xlsx format.", metavar="ANNOTATION"),
  make_option(c("-j", "--design"), default=NULL,
              help="Please write a DEseq2 design. Example: '~ treatment' or '~ treatment + cell_type' ", metavar="DESIGN"),
  make_option(c("-g", "--group"), default=NULL,
              help="Please specify group label for compairsion. Please choose one of column names in annotation table.", metavar="GROUP"),
  make_option(c("-t", "--treatment"), default=NULL,
              help="Please specify treatment name, please make sure it's same with [GROUP] column in the annotation table. DEG: log2(treatment/control)", metavar="TREATMENT"),
  make_option(c("-c", "--control"), default=NULL,
              help="Please specify control name, please make sure it's same with [GROUP] column in the annotation table. DEG: log2(treatment/control)", metavar="TREATMENT"),
  make_option(c("-O", "--output_path"), type="character", default=NULL,
              help="Please specify a output directory.", metavar="OUTPATH")
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Option validator
if (is.null(opt$counts)){
  print_help(opt_parser)
  stop("Please give a count table.", call.=FALSE)
}
if (is.null(opt$annotation)){
  print_help(opt_parser)
  stop("Please give a annotation table.", call.=FALSE)
}
if (is.null(opt$output_path)){
  print_help(opt_parser)
  stop("Please specify a output directory.", call.=FALSE)
}
if (is.null(opt$design)){
  print_help(opt_parser)
  stop("
        Please write a DEseq2 design. 

        Example: '~ treatment' or '~ treatment + cell_type' 
       
       ", call.=FALSE)
}
if (is.null(opt$group)){
  print_help(opt_parser)
  stop("
	Please specify a group label.

	Please choose one of column names in annotation table.	

       ", call.=FALSE)
}
if (is.null(opt$treatment) | is.null(opt$control)){
  print_help(opt_parser)
  stop("
        Please specify treatment and control name under [GROUP] column.

        Please make sure it's same with [GROUP] column in the annotation table!

	It will be calculate as log2(treatment/control).

       ", call.=FALSE)
}



# ----------------------------------------------------------------------------- #
print("Loading packages.")
## Load packages
pacman::p_load(
  rio,RColorBrewer,pheatmap,DESeq2,
  Hmisc,dplyr,ggplot2,DEGreport
)

# ----------------------------------------------------------------------------- #
# Loading data 
print("Loading data and annotation.")

data_format <- tools::file_ext(opt$counts)
raw_data <- import(opt$counts,format = data_format)
count_data <- raw_data[,-1]
rownames(count_data) <- raw_data[,1]

# Loading sample annotation!
anno_format <- tools::file_ext(opt$annotation)
anno.sample <- import(opt$annotation,format = anno_format)
rownames(anno.sample) <- anno.sample[,1]

# Check point
if (length(rownames(anno.sample))!=length(colnames(count_data))){
  stop("The number of samples in raw count table do not match the annotation table!", call.=FALSE)
}
checkpoint <- all(rownames(anno.sample)==colnames(count_data))
if (!checkpoint){
  stop("
	
	The order of ID in annotation and raw count table are different!
	
	Please edit the order of ID in annotation table.
	
	", call.=FALSE)
}

# ----------------------------------------------------------------------------- #
## Create DESeq2 Dataset object
print("Run DEseq2")

deseq_design <- as.formula(opt$design)
dds <- DESeqDataSetFromMatrix(
  countData = count_data, colData = anno.sample,
  design = deseq_design)

## Fit the model
dds <- DESeq(dds)

pdf(paste0(opt$output_path,"/","DEseq2_model.pdf"),width=8,height=6)
plotDispEsts(dds)
dev.off()

## counts normalize
normalized_counts <- data.frame(
  Gene=rownames(counts(dds, normalized=TRUE)),
  Normalized_count = counts(dds, normalized=TRUE)
)
export(normalized_counts,paste0(opt$output_path,"/","DEseq2_nomalized_count.csv"),format = "csv")

# ----------------------------------------------------------------------------- #

## =========================== Quality control =========================== ##
print("Calculate PCA stats.")
rld <- rlog(dds, blind=TRUE)

# calculate PCA 
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat),center = T)
# ----------------------------------------------------------------------------- #
# Each PC contribute how many % of variance? Scree plot
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

pdf(
  paste0(opt$output_path,"/","Scree Plot.pdf"),
  width = 6, height = 4
)

barplot(pca.var.per,
        names.arg = paste0("PC",seq(1,ncol(pca$x),1)),
        main="Scree Plot",
        ylim = c(0,as.integer(max(pca.var.per))+10),
        xlab="Principal Component",
        ylab="Percent Variation (%)",
        col = "lightcyan3",
        border = NA,
        cex.main = 1.5,
        cex.lab = 1,
        cex.axis = 1,
        font.lab = 1,
        font.axis = 1,
        yaxt="n",
        bty = "l"
)
# y-axis settings
axis(2,seq(0,as.integer(max(pca.var.per))+10,10),
     labels = T,
     lwd = 1,
     font = 1,
     cex.axis=1,
     las = 1
)

dev.off()

# ----------------------------------------------------------------------------- #

# ----------------------------------------------------------------------------- #
## Plot PCA using all genes

pdf(
  paste0(opt$output_path,"/","PCA all genes.pdf"),
  width = 8, height = 8
)

par(
  xpd=TRUE
)

## Calculate boundary
scalex <- (max(pca$x[,1]) - min(pca$x[,1]))*0.15
scaley <- (max(pca$x[,2]) - min(pca$x[,2]))*0.15
minPCAx <- as.integer(min(pca$x[,1])-scalex)
maxPCAx <- as.integer(max(pca$x[,1])+scalex)
minPCAy <- as.integer(min(pca$x[,2])-scaley)
maxPCAy <- as.integer(max(pca$x[,2])+scaley)

# plot PCA
plot(x=pca$x[,1], y=pca$x[,2],
     type = "p",
     main = "PCA",
     xlab = paste("PC1(",pca.var.per[1],"% variance)"),
     ylab = paste("PC2(",pca.var.per[2],"% variance)"),
     xlim = c(minPCAx,maxPCAx),
     ylim = c(minPCAy,maxPCAy),
     adj = 0.5,
     cex.main = 2.5,
     cex.lab = 1.5,
     cex.axis = 1.2,
     font.main = 2,
     font.lab = 2,
     ps = 2,
     pch = 16,
     cex = 1,
     col = "black",
     las = 1,
     bty = "l",
     tck = 1,
     tcl = -0.8,
     xaxs = "i",
     yaxs = "i"
)

cols <- brewer.pal(length(unique(anno.sample[,opt$group]))+1, "Set1")
points(x=c(pca$x[c(1,2),1],pca$x[c(1,2),1]),
       y=c(pca$x[c(1,2),2],pca$x[c(1,2),2]),
       cex = 1.5,
       pch = 16,
       col = cols[1]
)
points(x=c(pca$x[c(3,4),1],pca$x[c(3,4),1]),
       y=c(pca$x[c(3,4),2],pca$x[c(3,4),2]),
       cex = 1.5,
       pch = 16,
       col = cols[2]
)
xy <- par("usr")
#Figure legend
legend(
  x=xy[2]+xinch(-1.5),
  y=xy[3]+yinch(6.5),                               
  pch = 16,                                   
  col = cols,
  legend = unique(anno.sample[,opt$group]),
  bty = "b"
)

dev.off()

# ----------------------------------------------------------------------------- #
# Hack cluster tree!
# -----------------------------------------------------------------------------------------------------------------------------------------#
draw_dendrogram <- function(hc, gaps, horizontal = T) {
  # Define equal-length branches
  hc$height <- cumsum(rep(1/length(hc$height), length(hc$height)))
  h = hc$height/max(hc$height)/1.05
  m = hc$merge
  o = hc$order
  n = length(o)
  m[m > 0] = n + m[m > 0]
  m[m < 0] = abs(m[m < 0])
  dist = matrix(0, nrow = 2 * n - 1, ncol = 2, dimnames = list(NULL,
                                                               c("x", "y")))
  dist[1:n, 1] = 1/n/2 + (1/n) * (match(1:n, o) - 1)
  for (i in 1:nrow(m)) {
    dist[n + i, 1] = (dist[m[i, 1], 1] + dist[m[i, 2], 1])/2
    dist[n + i, 2] = h[i]
  }
  draw_connection = function(x1, x2, y1, y2, y) {
    res = list(x = c(x1, x1, x2, x2), y = c(y1, y, y, y2))
    return(res)
  }
  x = rep(NA, nrow(m) * 4)
  y = rep(NA, nrow(m) * 4)
  id = rep(1:nrow(m), rep(4, nrow(m)))
  for (i in 1:nrow(m)) {
    c = draw_connection(dist[m[i, 1], 1], dist[m[i, 2], 1],
                        dist[m[i, 1], 2], dist[m[i, 2], 2], h[i])
    k = (i - 1) * 4 + 1
    x[k:(k + 3)] = c$x
    y[k:(k + 3)] = c$y
  }
  x = pheatmap:::find_coordinates(n, gaps, x * n)$coord
  y = unit(y, "npc")
  if (!horizontal) {
    a = x
    x = unit(1, "npc") - y
    y = unit(1, "npc") - a
  }
  res = grid::polylineGrob(x = x, y = y, id = id)
  return(res)
}
# Replace the non-exported function `draw_dendrogram` in `pheatmap`:
assignInNamespace(x="draw_dendrogram", value=draw_dendrogram, ns="pheatmap")
# -----------------------------------------------------------------------------------------------------------------------------------------#
# ----------------------------------------------------------------------------- #
# plot Heatmap - connect to pheatmap script
rld_cor <- cor(rld_mat)    ## value in matrix --> color!

pdf(
  paste0(opt$output_path,"/","Correlation heatmap.pdf"),
  width = 10, height = 10
)

pheatmap(
  # Input matrix
  mat = rld_cor,
  
  # Style
  main = "Corelation heatmap",
  fontsize = 12,
  fontsize_row = 12,
  fontsize_col = 12,
  angle_col = 90,
  
  # Clustering
  cluster_rows = T,
  cluster_cols = T,
  clustering_distance_rows = "euclidean",
  #clustering_distance_cols = "euclidean",
  clustering_method =  "complete"
  # (options :"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".)
  # https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust
)

dev.off()

# ----------------------------------------------------------------------------- #
## ============================ Get DEGenes ============================ ##
## Define contrasts, extract results table, and shrink the log2 fold changes
contrast <- c(opt$group, opt$treatment, opt$control)
res_table_unshrunken <- results(
  dds, contrast=contrast, alpha = 0.05,independentFiltering = T)
res_table <- lfcShrink(
  dds, contrast=contrast, res=res_table_unshrunken,type = "normal")

## Result summary
summary(res_table)
res.df <- data.frame("Gene" = rownames(res_table),res_table)
export(res.df,paste0(opt$output_path,"/","DEG_Table.csv"),format = "csv")
# ----------------------------------------------------------------------------- #
## ============================= visualize ============================= ##
## MA plot
sigUP.df <- res.df[which(res.df$padj<0.05 & res.df$log2FoldChange>0),]
sigDOWN.df <- res.df[which(res.df$padj<0.05 & res.df$log2FoldChange<0),]

## XY
MApX <- log2(res.df$baseMean)
MApY <- res.df$log2FoldChange

## Exclude NA and Inf
filteredX <- MApX[!is.na(MApX) & is.finite(MApX)]
filteredY <- MApY[!is.na(MApY) & is.finite(MApY)]

## Calculate boundary
scaleMAx <- (max(filteredX) - min(filteredX))*0.15
scaleMAy <- (max(filteredY) - min(filteredY))*0.15
minMAx <- as.integer(min(filteredX)-scaleMAx)
maxMAx <- as.integer(max(filteredX)+scaleMAx)
minMAy <- as.integer(min(filteredY)-scaleMAy)
maxMAy <- as.integer(max(filteredY)+scaleMAy)

pdf(
  paste0(opt$output_path,"/","MAplot.pdf"),
  width = 10, height = 7
)

par(
  xpd=TRUE
)

plot(
  x = log2(res.df$baseMean),
  y = res.df$log2FoldChange,
  type = "p",
  main = "MA plot",
  xlab = "log2 mean expression",
  ylab = "log2 fold change",
  xlim = c(minMAx,maxMAx),
  ylim = c(minMAy,maxMAy),
  adj = 0.5,
  cex.main = 1.5,
  cex.lab = 1,
  cex.axis = 1.2,
  font.main = 2,
  font.lab = 2,
  ps = 1,
  pch = 16,
  cex = 0.8,
  col = "black",
  las = 1,
  bty = "l",
  tck = 1,
  tcl = -0.8,
  xaxs = "i",
  yaxs = "i"
)

points(
  x = log2(sigUP.df$baseMean),
  y = sigUP.df$log2FoldChange,
  cex = 0.8,
  pch = 16,
  col = "#B31B21"
)

points(
  x = log2(sigDOWN.df$baseMean),
  y = sigDOWN.df$log2FoldChange,
  cex = 0.8,
  pch = 16,
  col = "#1465AC"
)

xy <- par("usr")
#Figure legend
legend(
  x=xy[2]+xinch(-1.5),
  y=xy[3]+yinch(1),                               
  pch = 16,                                   
  col = c("#B31B21","#1465AC"),
  legend = c("UP regulated","DOWN regulated"),
  bty = "b"
)

dev.off()

##
## ============================= visualize ============================= ##
