###############################################################################
## Purpose: 
##    To design primers for targeted cDNA seq for fusion detection 
##    using anhcored multiplex PCR.
## 
## Reference:
##   Zheng Z, et al. Anchored multiplex PCR for targeted next-generation sequencing.
##   Nat Med. 2014 (http://www.nature.com/nm/journal/v20/n12/full/nm.3729.html)
##
## Contact     : Zongli Zheng, The Iafrate Lab, Massachusetts General Hospital
## Date created: 2012-04-13
##
## === Usage example: ===
##
## 1. Creat a panel gene list file with RefSeq transcript IDs "NM_###",
##      see the file 'example/lung.fusion.genelist.txt' for example.
##   
##    The NM numbers will be used to retreive template sequences from human
##    reference genome (GRCh37/hg19) for both sense and antisence sequences.
##
##    Alternatively, you can use your own template seqeunces (e.g. virus) under
##    the ~/project-AMP/$panel/seq/ folder (under development).
##
##  2. Run the command (see also 'example.command.txt')
##
##  > Rscript /path/to/AMP-primer-design/main.by.NM_cDNA.R \
##          assaytype=\"fusion"\
##          panel=\"lung.fusion\"  \
##          genelist=\"~/project-AMP/lung.fusion.genelist.txt\"  \ 
##          pjdir=\"~/project-AMP\"  \ 
##          depdir=\"~/dependency-AMP\"  \
##          ampdir=\"~/repo/AMP-primer-design\"  \ 
##          primer3path=\"~/primer3-2.3.5\"  \ 
##          blatdir=\"~/bin/x86_64\" \
##          keep_gfSvr=1 \
##          GSP1tag=\"GGATCTCGACGCTCTCCCT\"  \
##          GSP2tag=\"CCTCTCTATGGGCAGTCGGTGAT\"  \
##          ncpu=8  \
##          tempsize=90  \
##          NGSadaptors_and_humanRep=\"NGSadaptors_and_humanRep.fa\" 
##
##    where,
##     assaytype   - 'fusion' or 'mutation'. Fusion assay will retreive
##                   exonic sequence template for primer design, and 
##                   mutation assay intronic template.
##     panel       - the name of panel (e.g. \"lung.fusion\")
##     genelist    - name (with path) of the gene list, see 
##                     'example/lung.fusion.genelist.txt' for example.
##     pjdir       - is the project folder (required). A project/panel 
##                     folder e.g.  '~/project-AMP/lung.fusion' will
##                     be created by the pipeline
##     depdir      - path to dependency datasets
##     ampdir      - path to AMP primer design
##     primer3path - path to Primer3
##     blatpath    - path to BLAT
##     keep_gfSvr  - keep or free BALT gfSever in memory
##     GSP1tag     - the tag to be appended to 5' end of GSP1 primers.
##                     This tag do not participate in sequencing.
##     GSP2tag     - the tag to be appended to 5' end of GSP2 primers.
##                     For Illumina, this tag is Read2 Sequencing Primer.
##                     Default here is Ion Torrent (P23) sequence. This 
##                     tag allows for the same pirmers (hunreds to thousands)
##                     to be used for both Ion Torrent and Illumina platforms.
##                     (For Illumina Miseq, if use this GSP2 tag, in wet-lab:
##                        a. Add 3 ul of 100 μM of Illumina.custom.Index1.sequencing.primer
##                           to Miseq Reagent cartridge position 13 (Index Primer Mix)
##                        b. Add 3 ul of 100 μM of Illumina.custom.Read2.sequencing.primer
##                           to Miseq Reagent cartridge position 14 (Read 2 Primer Mix).
##                      See NGSadaptors.fa for the above primer seqeunces).
##     ncpu        - number of CPUs for multi-threading
##     tempsize    - size of template seqeunce to retreive from genome and
##                   to design primers on. Default 90 bp considers degraded
##                     RNA in FFPE samples
##
##     (Note: When running Rscript in Linux/Unix, remember to quote 
##        string values using \", and to add \ to separeate one command
##        into multiple lines, as shown in 'example.command.txt'.)
##
##
## === Computational steps of the pipeline ===
##
## 1. Generate template sequences (both sense and antisense)
##      and mask every dbSNPv137 locus as "N".
##
## 2. Use Primer3 to design candiate RIGHT primers.
##
##      - Up to 6 iterations of decreased design stringency are used to 
##      design candidate primers. Failed to design targets will be 
##      reported. creased stringency. 
##
##      - Primers are filtered against NGS adaptors including 
##      96 Illumina N5s, 12 N7s, 96 IonTorrent adaptors and 
##      human repetative sequences.
##
## 3. Map candidate primers against human genome
##    - Filter those mapped more than 5 (arbitrary) locations in the genome.   
##
## 4. Pairing candidates for GSP1 and GSP2 based on pair penalty scores.
##
## 
## === An example - a lung fusion panel ===
##
## Please see the example folder for how to create a list of genes, and if
##   needed, exons and sense/antisense primers, as well as final result
##   output on primer sequences, primer bed file etc.
##
###############################################################################

# get initial values
args=(commandArgs(TRUE))
for(i in 1:length(args)){eval(parse(text=(args[[i]])))}

## create working dir for the panel
paneldir = paste(pjdir, panel, sep='/')
dir.create(paneldir, showWarnings = FALSE, recursive = TRUE)

setwd(paneldir) 

## get gene list
gene.list = read.table(genelist, sep='\t', header=T, fill=T)

## remove blank space, heading 0s
gene.list$NM = gsub(' ', '', gene.list$NM)
gene.list$exon = tolower(gsub(' ', '', gene.list$exon))
gene.list$exon = gsub('^0+', '', gene.list$exon)

if (assaytype == 'fusion'){
	gene.list$sense = tolower(gsub(' ', '', gene.list$sense))
} else {
	## for compatability with mutation assay, assign sense 'both'
	gene.list$sense = 'both'
}

panel.NM = sort(unique(gene.list$NM))
writeLines(panel.NM, 'panel.NM')
system(paste("join panel.NM ", depdir, "/hg19.RefGene.NM > target.refseq", sep=''))

######################################################################
# step 1. retreive sequence from genome 
source(paste(ampdir, '/step.1.retreive.sequence.R', sep='')) 
	

######################################################################
# step 2. call primer3 - several iteration rounds
source(paste(ampdir, '/step.2.call.primer3.R', sep='')) 


######################################################################
# setp 3. BLAT candidate primers against genome
	t = 11
	s = 4
	n.srv = ceiling(ncpu/4)
	ClnPsrv = 3
	n.fa=n.srv*ClnPsrv
	repMatch=1024
	start.port=9000
	faDir = 'split'
	depdir
	pjdir
	db = 'hg19.2bit'
	minIden = 93
	minScore = 12
	maxIntron = 700001
	NAtype = 'dna'
source(paste(ampdir, '/step.3.blat.candidate.primer.R', sep='')) 

# organize balt results
system(paste('bash ', ampdir, '/scripts/afterBlat.sh', sep=''))


######################################################################
# step 4. pair GSP1-GPS2
source(paste(ampdir, '/step.4.pairing.GSPs.R', sep='')) 

	## top candidates
	system('head -n1 ranked.all.candidate.pairs > ranked.1st.pairs')
	system("grep ' 1$' ranked.all.candidate.pairs >> ranked.1st.pairs")


######################################################################
# step 5. check uniqueness of 12 bases at 3' of all primers. 
#      - If not unique, candidate pairs from 'ranked.all.candidate.pairs'
#	    will be retreived, ranked 2/3/.., until all tail 12 bases
#           are unique.
#
#      - If all unique, proceed.
#
#      - Save final tail 12 bases of all primers for future use.
#          e.g, to check uniqueness when adding new primers 
#	   to an existing panel.
#
pairs1 = read.table('ranked.1st.pairs', header=T, stringsAsFactors=F)
pairs1$target = sub(":.*", "", pairs1$r1.qName)
source(paste(ampdir, '/step.5.check.uniqueness.R', sep='')) 

writeLines(all.t12, paste(panel, '.tail.12bases.'
	  , format(Sys.time(), "%Y.%b.%d"), sep=''))
	   

######################################################################
# Lastly,
#  - add GSP tags
#  - output 
#  - save primer.bed
######################################################################
pairs.keep$GSP1 = paste(GSP1tag, pairs.keep$r1.seq,sep='')
pairs.keep$GSP2 = paste(GSP2tag, pairs.keep$r2.seq,sep='')
pairs.keep$gene = sapply(strsplit(pairs.keep$target, '_'), "[[", 1)
pairs.keep$exon = sapply(strsplit(pairs.keep$target, '_'), "[[", 2)
pairs.keep$exonSize = sapply(strsplit(pairs.keep$target, '_'), "[[", 3)
if (assaytype == 'fusion'){
	pairs.keep$sense = sapply(strsplit(pairs.keep$target, '_'), "[[", 4)
} else {
	pairs.keep$sense = pairs.keep$exonSize
}

gene.NM = ref[,c('name2','name')]
names(gene.NM) = c('gene', 'NM')
final = merge(pairs.keep, gene.NM, by='gene')
final$gsp1.name = paste(final$gene, '_ex', final$exon, '_', final$sense
			  , '.1 (', final$NM, ')', sep='')
final$gsp2.name = paste(final$gene, '_ex', final$exon, '_', final$sense
			  , '.2 (', final$NM, ')', sep='')

# output all exons
	gsp1 = final[, c('gsp1.name', 'GSP1')]
	gsp2 = final[, c('gsp2.name', 'GSP2')]

	# primer bed - only GSP2 is sequenced/relevant
	primer.bed = final[, c('r1.chr', 'r2.tStart', 'r2.tEnd', 'gsp2.name')]

	write.csv(gsp1, paste(panel, '_all.gsp1.csv', sep=''), quote=F, row.names=F)
	write.csv(gsp2, paste(panel, '_all.gsp2.csv', sep=''), quote=F, row.names=F)
	write.table(primer.bed, paste(panel, '_all.gsp2.primer.bed', sep=''), sep='\t'
		    , quote=F, row.names=F, col.names=F)

# selected exon/sense
	if (assaytype == 'fusion'){
		final$nm.sense = paste(final$NM, tolower(final$sense), sep='_')
		gene.list$nm.sense = paste(gene.list$NM, gene.list$sense, sep='_')
		nm.sense = unique(gene.list$nm.sense)
		final$sense.select = 0
		final$sense.select[final$nm.sense %in% nm.sense] = 1
	}

	final$nm.exon = paste(final$NM, as.numeric(final$exon), sep='_')
	gene.list$nm.exon = paste(gene.list$NM, gene.list$exon, sep='_')
	nm.exon = unique(gene.list$nm.exon)
	final$exon.select = 0
	final$exon.select[final$nm.exon %in% nm.exon] = 1
	nm.all.exons = gene.list$NM[gene.list$exon=='all']
	final$exon.select[final$NM %in% nm.all.exons] = 1

	final$select = 0
	if (assaytype == 'fusion'){
		final$select[final$exon.select ==1 & final$sense.select ==1] = 1
	} else {
		final$select[final$exon.select ==1] = 1
	}


	## backup final data
	 write.table(final, paste(panel, '_intermediate.data.txt', sep=''), sep='\t'
		     , quote=F, row.names=F)
	

	final.s = subset(final, select==1)

	##############################
	## final output
	##############################
	GSP1 = final.s[, c('gsp1.name', 'GSP1')]
	GSP2 = final.s[, c('gsp2.name', 'GSP2')]

	# primer bed - only GSP2 is sequenced/relevant
	primer.bed.s = final.s[, c('r1.chr', 'r2.tStart', 'r2.tEnd', 'gsp2.name')]

	write.csv(GSP1, paste(panel, '_GSP1.csv', sep=''), quote=F, row.names=F)
	write.csv(GSP2, paste(panel, '_GSP2.csv', sep=''), quote=F, row.names=F)
	write.table(primer.bed.s, paste(panel, '_GSP2.primer.bed', sep=''), sep='\t'
		    , quote=F, row.names=F, col.names=F)


## cleanup
	system('mv unpaired.targets  targets.failed.pairing.GSPs.txt')
	if (file.exists('missed.seq_seq_primer3.exome.setting.3')){
		system('mv missed.seq_seq_primer3.exome.setting.3  targets.failed.primer3.design.txt')
	}
	system('rm -rf split seq seq.noMask out all.psl.matt.sorted primer.pos.tm.sorted primer.psl missed.seq_* splitFa.Rout panel.NM primer.candidates.0 primer3.exome.setting*')

## END

