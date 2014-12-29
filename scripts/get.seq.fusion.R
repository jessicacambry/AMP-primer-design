
## This script is to be called by 'step.1.retreive.sequence.R' to:
##   - retreive exonic seqeuence templates from genome (hg19/b37) 
##       for detection of 5' and 3' fusions.

# === Requires:
#  gene - the target
#  tempsize - size of template sequences to retreive
#  depdir - dependency data path

# === Should have existed:
#  seq/
#  seq.noMask/
#  target.refseq

## usage e.g.: 
##   Rscript get3n5AMPseq.R gene tempsize depdir


library(BSgenome.Hsapiens.UCSC.hg19)
library(SNPlocs.Hsapiens.dbSNP.20111119)
ref = read.table('target.refseq',header=F,stringsAsFactors=F)
names(ref)=c('name',    'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames')

gene = commandArgs(TRUE)[1]
tempsize = as.numeric(commandArgs(TRUE)[2])
depdir = commandArgs(TRUE)[3]

system(paste("touch _running_", gene, sep=''))

if (gene %in% ref$name2){

   print(gene)
   # make exon label strings 001, 002...
   iii = substr(1001:1999, 2, 5)

   geneTab = ref[ref$name2 == gene,]
   chrom = sub('chr','',geneTab$chrom)
   snp.chr = read.table(paste(depdir, '/dbsnp/snp.', chrom,sep=''), stringsAsFactors=F)

   start = as.numeric(unlist(strsplit(as.character(geneTab$exonStarts),','))) -1
   end = as.numeric(unlist(strsplit(as.character(geneTab$exonEnds),','))) +1
   strand = geneTab$strand
   n.exon = geneTab$exonCount
   size = end - start
  
   domain.start = 1
   domain.end = n.exon

   ## +++++++++++++++++++++++++++++++++++++++++++++++++
   if (strand=='+'){
	   ex1 = 1
	   ex2 = n.exon
	   ex3 = 1
	   ex4 = n.exon
   }

   if (strand=='-'){
	   ex1 = n.exon
	   ex2 = 1
	   ex3 = n.exon 
	   ex4 = 1
   }
   ## --------------------------------------------------

   chrom = paste('chr', chrom, sep='')
   ## 12121212121212121212121212121212121212
   for (ex in ex1:ex2){
	   start.ex = start[ex]
	   size.ex = size[ex]
	   ex.number = ifelse(strand=='+', ex, n.exon-ex+1)
	   getSeq = data.frame('getSeq'=getSeq(Hsapiens, chrom, start.ex+2+5
	      , start.ex+min(size.ex-1, tempsize)-1, as.character=T, strand='+')
	      , 'size'=size.ex, 'exon'=ex.number)

	   getSeq$seq.noMask = getSeq$getSeq
		   ## mask dbsnp
		   snp12 = subset(snp.chr, V2 >= start.ex+7 & V2 <= start.ex+min(size.ex, tempsize))
		   snp12n = length(snp12$V1)
		   seq12 = getSeq$getSeq
		   if (snp12n>0){
			   for (i in 1:snp12n){
			   pos = snp12$V2[i] - (start.ex+7) +1
			   Ns = nchar(snp12$V4[i])
			   l.seq = substr(seq12, 1, pos -1)
			   r.seq = substr(seq12, pos + Ns, min(size.ex,tempsize)+1)
			   seq12 = paste(l.seq, paste(rep('N', Ns), collapse=''), r.seq, sep='')
			   }
		   getSeq$getSeq = seq12
		   }
		   
		   ## output
		   race = ifelse(strand=='+', '_AntiSense', '_Sense')
		   getSeq$getSeq = paste('SEQUENCE_ID=', gene,'_',iii[ex.number],'_',getSeq$size, getSeq$strand, race, '\nSEQUENCE_TEMPLATE=',getSeq$getSeq,'\n=',sep='')
		   getSeq$seq.noMask = paste('SEQUENCE_ID=', gene,'_',iii[ex.number],'_',getSeq$size, getSeq$strand, race, '\nSEQUENCE_TEMPLATE=',getSeq$seq.noMask,'\n=',sep='')
		   write.table(getSeq$getSeq, file=paste('seq/', gene,'_',iii[ex.number],'_',getSeq$size, getSeq$strand, race, sep=''),row.names=F,quote=F,col.names=F)
		   write.table(getSeq$seq.noMask, file=paste('seq.noMask/', gene,'_',iii[ex.number],'_',getSeq$size, getSeq$strand, race, sep=''),row.names=F,quote=F,col.names=F)
	   }

   ## 3434343434343434343434343434343434343434
   for (ex in ex3:ex4){
	   end.ex = end[ex]
	   size.ex = size[ex]
	   ex.number = ifelse(strand=='+', ex, n.exon-ex+1)
	   getSeq = data.frame('getSeq'=getSeq(Hsapiens, chrom, end.ex-min(size.ex-1,tempsize)+1
	      , end.ex -1-5, as.character=T,strand='-')
	      , 'size'=size.ex, 'exon'=ex.number)

	   getSeq$seq.noMask = getSeq$getSeq
		   ## mask dbsnp
		   snp34 = subset(snp.chr, V2 >= (end.ex-min(size.ex-1,tempsize)+1) & V2 <= end.ex-6)
		   snp34n = length(snp34$V1)
		   seq34 = getSeq$getSeq
		   if (snp34n>0){
			   for (i in 1:snp34n){
			   pos = end.ex-1-5 - snp34$V2[i] 
			   Ns = nchar(snp34$V4[i])
			   l.seq = substr(seq34, 1, pos)
			   r.seq = substr(seq34, pos + Ns +1, min(size.ex-1, tempsize)+2)
			   seq34 = paste(l.seq, paste(rep('N',Ns),collapse=''), r.seq, sep='')
			   }
		   getSeq$getSeq = seq34
		   }
		   
		   ## output
		   race = ifelse(strand=='+', '_Sense', '_AntiSense')
		   getSeq$getSeq = paste('SEQUENCE_ID=', gene,'_',iii[ex.number],'_',getSeq$size, getSeq$strand, race, '\nSEQUENCE_TEMPLATE=',getSeq$getSeq,'\n=',sep='')
		   getSeq$seq.noMask = paste('SEQUENCE_ID=', gene,'_',iii[ex.number],'_',getSeq$size, getSeq$strand, race, '\nSEQUENCE_TEMPLATE=',getSeq$seq.noMask,'\n=',sep='')
		   write.table(getSeq$getSeq, file=paste('seq/', gene,'_',iii[ex.number],'_',getSeq$size, getSeq$strand, race, sep=''),row.names=F,quote=F,col.names=F)
		   write.table(getSeq$seq.noMask, file=paste('seq.noMask/', gene,'_',iii[ex.number],'_',getSeq$size, getSeq$strand, race, sep=''),row.names=F,quote=F,col.names=F)
	   }
} else {
	   cat(gene, " not found, skipped...\n") 
	   cat(gene, file = "GENE_SYMBOL_NOT_FOUND.txt", sep='\n', append=T) 
}

system(paste("rm _running_", gene, sep=''))
