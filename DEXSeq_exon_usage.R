### DEXSeq Differential Exon Usage Analysis
# initial steps of analysis are done using two Python scripts provided with DEXSeq
#
# pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )
# list.files(pythonScriptsDir)

# system.file( "python_scripts", package="DEXSeq", mustWork=TRUE )

###  Preparing the annotation
# ubuntu command line shell
# Make sure that your current working directory contains the GTF file 
# python path_to_script GTF_file_name output_file.gff
# ex. python /mnt/c/Users/alman/Documents/R/win-library/3.6/DEXSeq/python_scripts/dexseq_prepare_annotation.py Drosophila_melanogaster.BDGP5.72.gtf Dmel.BDGP5.25.62.DEXSeq.chr.gff

###  Counting reads per exon
# for each sam/bam file
# python /path/to/dexseq_count.py output_file.gff untreated1.sam untreated1fb.txt 
# -f bam
# -s reverse

### Reading the data in to R 
library("DEXSeq")
library("reshape2")
library("ggplot2")
library("stringr")

setwd("C:/Users/alman/Desktop/RNAseq")
countFiles = list.files( pattern=".txt$", full.names=TRUE)
countFiles <- countFiles[c(1:3,7:9)]
flattenedFile = list.files( pattern="gff$", full.names=TRUE)

sampleTable = data.frame(
  row.names = c( "D24R1", "D24R2", "D24R3", 
                 "M24R1", "M24R2", "M24R3" ),
  condition = c("DMSO", "DMSO", "DMSO",  
                "MKC", "MKC", "MKC" ),
  libType = c( "single-end", "single-end", "single-end", 
               "single-end", "single-end", "single-end" ))

#construct an DEXSeqDataSet object from this data
dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )

# first 6 columns = number of reads mapping to exonic regions, last 6 = sum of the counts for rest of the exons from the same gene on each sample.

### Normalization (same method as in DEseq2)
dxd = estimateSizeFactors( dxd )

# Dispersion estimation (estimate strength of the noise)
dxd = estimateDispersions( dxd )
plotDispEsts( dxd )

### Testing for differential exon usage
# fits a generalized linear model with the formula ~sample + exon + condition:exon and compare it to the smaller model (the null model) ~ sample + exon.

dxd = testForDEU( dxd )

# estimate relative exon usage fold changes
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")

### Results summary
dxr1 = DEXSeqResults( dxd )
mcols(dxr1)$description
resultados <- as.data.frame(dxr1)

### VISUALIZATION
# MA plot
plotMA( dxr1, cex=0.8 )
table ( dxr1$padj < 0.1 )
hits <- subset(resultados, resultados$padj < 0.1)
hits_table <- as.data.frame(table(hits$groupID))

# plot exon usage
plotDEXSeq( dxr1, "ENSG00000026508", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

### Save summary results in HTML bundle
DEXSeqHTML(dxr1, FDR = 0.1)

#plot frecuencies pie
exon_frequency_pie <- function(nombre_gen) {
  nombre_gen1  <- subset(resultados, resultados$groupID == nombre_gen)
  m1_sums <- na.omit(nombre_gen1$DMSO)/log2(na.omit(nombre_gen1$genomicData.width))
  m2_sums <- na.omit(nombre_gen1$MKC)/log2(na.omit(nombre_gen1$genomicData.width))
  par(mfrow = c(1,2))
  pie(m1_sums, main="DMSO")
  pie(m2_sums, main="MKC")
}  

exon_frequency_pie("ENSG00000136942")

### GET SEQUENCES OF THE EXONS
### add another column with the nucleotide seqs
## might want to look at range around it, not just there
library(BSgenome.Hsapiens.NCBI.GRCh38)

chromO <- hits$genomicData.seqnames
starT <- hits$genomicData.start
enD <- hits$genomicData.end

sequence_exon <- function(x) {
  as.character(getSeq(Hsapiens, chromO[x], start = starT[x], end = enD[x]))
}

sequence_exon_flank <- function(x) {
  as.character(getSeq(Hsapiens, chromO[x], start = starT[x]-500, end = enD[x]+500))
}

for(i in 1:nrow(hits)) {
  hits$sequence[i] <- sequence_exon(i)
}

for(i in 1:nrow(hits)) {
  hits$sequence_flanks[i] <- sequence_exon_flank(i)
}

### scan for consensus cleavage sequence in exon seqs
consensus_seq <- as.character("CTGCAG")

for(i in 1:nrow(hits)) {
  hits$cleavage[i] <- str_detect(hits$sequence[i], consensus_seq)
}

for(i in 1:nrow(hits)) {
  hits$cleavage_flank[i] <- str_detect(hits$sequence_flanks[i], consensus_seq)
}

############################################
### Screen ALL exons for a given gene
### TAKE A LOOK AT EACH GENE AT A TIME
library("biomaRt")

listEnsembl(GRCh=37)
listEnsembl(version=91)

ensembl = useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

# LOOK FOR CONSENSUS CLEAVAGE SEQUENCE trough all exons of one gene
CD44_exons <- getSequence( id="CD44" ,type='hgnc_symbol',seqType = 'gene_exon', mart = ensembl)
# CD44 <- getSequence(chromosome = 4, start = 50000-3000, end = 50000+3000,type='hgnc_symbol',seqType = 'gene_exon_intron', mart = GENES)
CD44_exons$cleavage <- "NA"

for(i in 1:nrow(CD44_exons)) {
  CD44_exons$cleavage[i] <- str_detect(CD44_exons$gene_exon[i], consensus_seq)
}

for(i in 1:nrow(CD44_exons)) {
  CD44_exons$lenght[i] <- nchar(CD44_exons$gene_exon[i])
}

for(i in 1:nrow(CD44_exons)) {
  CD44_exons$residues[i] <- nchar(CD44_exons$gene_exon[i])/3
}

###
CD44_exons <- subset(resultados, resultados$groupID == "ENSG00000026508")

for(i in 1:nrow(CD44_exons)) {
  CD44_exons$sequence[i] <- sequence_exon(i)
}

for(i in 1:nrow(CD44_exons)) {
  CD44_exons$cleavage[i] <- str_detect(CD44_exons$sequence[i], consensus_seq)
}

### annotate gene symbol

probando <- function(x) {
getSequence(chromosome = hits_seq_cleav_annot$genomicData.seqnames[x], start = hits_seq_cleav_annot$genomicData.start[i], end = hits_seq_cleav_annot$genomicData.end[i],type='hgnc_symbol',seqType = 'gene_exon_intron', mart = ensembl)[2]
}

for(i in 1:nrow(hits_seq_cleav_annot)) {
  hits_seq_cleav_annot$hgnc_symbol[i] <- probando(i)
}
