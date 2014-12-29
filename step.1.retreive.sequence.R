##
## To be called by main.R to generate template sequences
##   and mask every dbSNPv137 locus as "N" before primer design. 
##   Different sequences will be retreived from genome (hg19/CRCh37) 
##   for different assay types:
##
##   - For 'fusion' assay, exonic sequences (both sense and antisense)
##        will be retreived.
##
##   - For 'mutation' assay, which detects coding seqeunce mutation 
##        (including SNV, Indel and CNV):
##
##   	  -- For most exons (smaller than predifined 'subExonSize'),
##             two template sequences will be retreived.
##
##            1) a sequence of length 'tempsize' at 3' of the exon on
##               the plus (forward) strand, with a 'leadsize' distance
##               to the exon.
##            2) a sequence of length 'tempsize' at 5' of the exon on
##               the minus (reverse) strand, with a 'leadsize' distance
##               to the exon.
##            One primer set for each of the strands will be designed
##              for bi-template sequencing coverage.
##
##        -- For eoxons larger than 'subExonSize', tiling of template 
##             sequences will be retreived for tiled coverage across
##             the large exons (also bi-template).
## 
## requires:
##   - assaytype
##   - target.refseq
##   - subExonSize (for mutation assay)

ref = read.table('target.refseq', header=F, stringsAsFactors=F)
names(ref)=c('name',    'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames')
# head(ref, 3)

## retreive sequences by refseq NM (by gene names or by chr.pos not supported yet)
genes = ref$name2

# clear files if exist
system('rm -rf seq/ seq.noMask/', ignore.stderr=F)
system('mkdir seq seq.noMask')

# multi-threading
nrun=length(list.files('./', '_running_'))
if (nrun>0){system('rm _running_*')}
nrun = 1


for (gene in genes){
        while (nrun >= ncpu){
                Sys.sleep(5)
                nrun=length(list.files('./', '_running_'))
        }

	if (assaytype == 'fusion'){
		# fusion assay
		system(paste("Rscript --vanilla ", ampdir, "/scripts/get.seq.fusion.R ",
			     gene, " ", tempsize, " ", depdir
			     , " &", sep=''))
	} else if (assaytype == 'mutation'){
		## mutation assay
		system(paste("Rscript --vanilla ", ampdir, "/scripts/get.seq.mutation.R ",
			     gene, " ", tempsize, " ", depdir, " ", subExonSize, " "
			     , leadsize, " &", sep=''))
	} else {stop ("assaytype should be 'fusion' or 'mutation'")}
}

        # wait until all genes are done
        Sys.sleep(10)
        while (nrun >0){
                Sys.sleep(5)
                nrun=length(list.files('./', '_running_'))
        }

## End
