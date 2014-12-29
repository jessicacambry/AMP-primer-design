
## This script is to be called by 'step.1.retreive.sequence.R' to:
##   - retreive sequences from genome (hg19/GRCh37) for detection of 
##     coding seqeunce mutation (including SNV, Indel and CNV).
##   


# === Requires:
#  gene - the target
#  tempsize - size of template sequences to retreive
#  depdir - dependency data path
#  subExsize - for exons larger than which, design titling primers
#  leassize - the distance between retreive sequence template 
#               (normally intronic) and target (normally exons)
#               

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
subExonSize = as.numeric(commandArgs(TRUE)[4]) 
leadsize = as.numeric(commandArgs(TRUE)[5])

system(paste("touch _running_", gene, sep=''))

if (gene %in% ref$name2){

   print(gene)
   # make exon label strings 001, 002...
   iii = substr(1001:1999, 2, 5)

   geneTab = ref[ref$name2==gene,]
   chrom = geneTab$chrom
   strand = geneTab$strand
   start0 = as.numeric(unlist(strsplit(as.character(geneTab$exonStarts),','))) -1
   end0 = as.numeric(unlist(strsplit(as.character(geneTab$exonEnds),','))) +1
   size0 = end0 - start0
   n.exon0 = length(end0)

   ## read dbsnp
   chrn = toupper(sub('chr', '', chrom))
   snp.chr = read.table(paste(depdir, '/dbsnp/snp.', chrn, sep=''), stringsAsFactors=F)
   
   # default to ignore utr sequnece and design primer for coding sequences only.
   if (!exists('utr')){utr = 0}
   if (utr=='0'){
        start = ifelse(start0 >= geneTab$cdsStart, start0, geneTab$cdsStart)
        end = ifelse(end0 <= geneTab$cdsEnd, end0, geneTab$cdsEnd)
   }else{
        start = start0
        end = end0
   }
   size = end - start
   exs = which(size >0)
   n.exon = length(exs)

  for (ex in exs){
           ex.sub.n = ceiling(size[ex] / subExonSize)
           ex.sub.size = ifelse(ex.sub.n == 1, 0, floor(size[ex]/ex.sub.n) - 1)
           ex.i = ifelse(strand == '+', ex, n.exon0 - ex + 1)

       for (ex.sub in ex.sub.n:1){
	   ########################################
	   ### 3' of exon on plus (fwd) strand
	   ########################################
	   end.ex.sub = end[ex] - ex.sub.size * (ex.sub - 1) 
			+ leadsize # at least 'leadsize' bp in intron
	   end.ex.sub.tempsize = end.ex.sub + tempsize
	   end.ex.sub2 = ifelse(ex == rev(exs)[1], end.ex.sub.tempsize,
				min(end.ex.sub.tempsize, start[ex + 1] 
				  - leadsize)
				)
	   end.ex.sub.TplSize = end.ex.sub2 - end.ex.sub


    	if (end.ex.sub.TplSize > 30){
	   ex.3f = data.frame('getSeq'=getSeq(Hsapiens, chrom
			, end.ex.sub, end.ex.sub2, as.character=T
			, strand='+')
			, 'size' = size[ex]
			, 'sub' = letters[ex.sub]
			, 'strand'='+')

	   ## before mask dbsnp
	   ex.3f$seq = paste('SEQUENCE_ID=', gene, '_', iii[ex.i]
				, '_', ex.3f$size, ex.3f$sub
				, ex.3f$strand
			, '\nSEQUENCE_TEMPLATE=', ex.3f$getSeq
			, '\n=', sep='')
	   write.table(ex.3f$seq, file=paste('seq.noMask/', gene
			, '_', iii[ex.i], '_', ex.3f$size, ex.3f$sub
			, ex.3f$strand, sep='')
			, row.names=F, quote=F, col.names=F)

	   ## mask dbsnp
	   snp = subset(snp.chr, V2 >= end.ex.sub 
				& V2 <= end.ex.sub + tempsize)
	   snpn = length(snp$V1)
	   seq = ex.3f$getSeq
	   if (snpn>0){
		   for (i in 1:snpn){
		   pos = snp$V2[i] - end.ex.sub +1
		   Ns = nchar(snp$V4[i])
		   l.seq = substr(seq, 1, pos -1)
		   r.seq = substr(seq, pos + Ns, end.ex.sub +
				 tempsize + 1)
		   seq = paste(l.seq, paste(rep('N', Ns), collapse='')
				, r.seq, sep='')

			}
	   ex.3f$getSeq = seq
	   }

	   ex.3f$seq = paste('SEQUENCE_ID=', gene, '_', iii[ex.i]
			, '_', ex.3f$size, ex.3f$sub, ex.3f$strand
			, '\nSEQUENCE_TEMPLATE=', ex.3f$getSeq
			, '\n=', sep='')
	   write.table(ex.3f$seq, file = paste('seq/', gene, '_'
			, iii[ex.i], '_', ex.3f$size, ex.3f$sub
			, ex.3f$strand, sep='')
			, row.names=F, quote=F, col.names=F)

    	} else {
		cat(c(gene, iii[ex.i], strand, ex, ex.sub, end.ex.sub.TplSize
			, '\n'), file = 'smallTplSize', append=T)
               }

		
	   ########################################
      	   # 5' of exon on the minus (revers) strand
	   ########################################
	   start.ex.sub = start[ex] + ex.sub.size * (ex.sub - 1) 
				- leadsize 
	   start.ex.sub.tempsize = start.ex.sub - tempsize
	   start.ex.sub2 = ifelse(ex == exs[1], start.ex.sub.tempsize,
				  max(start.ex.sub.tempsize, end[ex - 1]
					 + leadsize)
				  )
	   start.ex.sub.TplSize = start.ex.sub - start.ex.sub2

	   if (start.ex.sub.TplSize > 30){

	   ex.5r = data.frame('getSeq' = getSeq(Hsapiens, chrom
			, start.ex.sub2, start.ex.sub, as.character=T
			, strand = '-')
			, 'size' = size[ex]
			, 'sub' = letters[ex.sub]
			, 'strand' = '-')

	   ## before mask dbsnp
	   ex.5r$seq = paste('SEQUENCE_ID=', gene, '_', iii[ex.i]
			, '_', ex.5r$size, ex.5r$sub, ex.5r$strand
			, '\nSEQUENCE_TEMPLATE=', ex.5r$getSeq
			, '\n=', sep='')
	   write.table(ex.5r$seq, file = paste('seq.noMask/', gene
			, '_', iii[ex.i], '_', ex.5r$size, ex.5r$sub
			, ex.5r$strand, sep='')
			, row.names = F, quote = F, col.names=F)

	   ## mask dbsnp
	   snp5 = subset(snp.chr, V2 >= (start.ex.sub2) 
				& V2 <= start.ex.sub)
	   snp5n = length(snp5$V1)
	   seq5 = ex.5r$getSeq
	   if (snp5n>0){
		   for (i in 1:snp5n){
		   pos = start.ex.sub - snp5$V2[i] +1
		   Ns = nchar(snp5$V4[i])
		   l.seq = substr(seq5, 1, pos - Ns)
		   r.seq = substr(seq5, pos + 1, start.ex.sub.TplSize)
		   seq5 = paste(l.seq, paste(rep('N', Ns), collapse='')
				, r.seq, sep='')
		   }
	   ex.5r$getSeq = seq5
	   }
	   ex.5r$seq = paste('SEQUENCE_ID=', gene, '_', iii[ex.i], '_'
			, ex.5r$size, ex.5r$sub, ex.5r$strand
			, '\nSEQUENCE_TEMPLATE=', ex.5r$getSeq
			, '\n=',sep='')
	   write.table(ex.5r$seq, file = paste('seq/', gene, '_', iii[ex.i]
			, '_', ex.5r$size, ex.5r$sub, ex.5r$strand, sep = '')
			, row.names = F, quote = F, col.names = F)

	   }else{
	   	cat(c(gene, iii[ex.i], strand, ex, ex.sub, end.ex.sub.TplSize
			, '\n'), file = 'smallTplSize', append=T)
	   }
}

######################################
### write exon bed + (leadsize - 2)
######################################

        bed0 = paste(chrom, start[ex] - leadsize + 1, end[ex] + leadsize - 1
			, gene, ex.i, size[ex], ex.sub.n, sep = '\t')
        write.table(bed0, 'target.bed0', append = T, quote = F, col.names = F
			, row.names = F)
   }
}


system(paste("rm _running_", gene, sep=''))
