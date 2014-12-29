##
## This program is to be called by pairing.GSPs.R for multi-threading 
##   of pairing GSPs, by doing a subset of targets using one cpu
##

cpu = commandArgs(TRUE)[1]
ampdir = commandArgs(TRUE)[2]
setwd(cpu)

# get primer.psl
system("join targets ../primer.psl.srt > primer.psl.cpu")

targets0 = readLines('targets')
psl = read.table('primer.psl.cpu', stringsAsFactors=F)
names(psl) = c('tar', 'candName', 'pos', 'primerLen', 'TM', 'MIS'
	       , 'candSeq', 'matches', 'misMatches', 'repMatches', 'nCount'
	       , 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert'
	       , 'strand', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart'
	       , 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts'
	       , 'posLen')
psl$qName = paste(psl$tar, psl$candName, sep=':')
psl$head = psl$pos - psl$primerLen + 1

## select goodHits -- the lower number of hits to genome, the better
hitTab = data.frame(table(psl$qName))
names(hitTab) = c('cand', 'hitN')
hitTab$tar = sub(':.*', '', hitTab$cand)
# head(hitTab)

###########################
## select pairs from candidates with no more than 'hit' number of hits
###########################

## the function
# requires: psl hitTab
pairing.max.hit = function(tars, hitn, pslD, hitTabD) {
	goodHits = hitTabD$cand[hitTabD$hitN <= hitn]
	pslGood = subset(pslD, qName %in% goodHits)
	pslGood.tar = subset(pslGood, tar %in% tars)
	nrow(pslD); nrow(pslGood.tar)

	if (nrow(pslGood.tar) >= 2){
		good.hit.tars = unique(pslGood.tar$tar)

		## pairing based on pair-panelty, within a target
		for (good.hit.tar in good.hit.tars){
			good.hit.tar.psl = subset(pslGood.tar, tar == good.hit.tar)
			if (nrow(good.hit.tar.psl) >= 2){
				dput(good.hit.tar.psl, '__good.hit.tar.psl')
				source(paste(ampdir, "/scripts/pair.penalty.R", sep='')) 
				system('rm __good.hit.tar.psl')
			}
		}
		system(paste("cat __paired.cands.ranked_* > _paired.ranked.", hitn
			     , sep=''), ignore.stderr=T)
		system("rm __paired*", ignore.stderr=T)
		system(paste("cut -d ' ' -f10 _paired.ranked.", hitn, " | grep -v qName | sed 's/:.*//' | sort -u > _paired.targets.", hitn, sep=''))

		paired.targets = readLines(paste('_paired.targets.', hitn, sep=''))
		missed.targets = targets[!(targets %in% paired.targets)]
		} else {
			missed.targets=tars
	}
	writeLines(missed.targets, paste('_missed.targets.', hitn, sep=''))
}

# initialization
system('rm _*', ignore.stderr=T)

## iterations 
# round 1: hit=1
hit.min = 1
hit.max = 5
n.missed.tar=1
targets = targets0
for (hiti in hit.min:hit.max){
	if (n.missed.tar >0){

		## the pairing function
		cat('iteration: hit =', hiti, '\n')
		pairing.max.hit(tars=targets, hitn=hiti, pslD=psl, hitTabD=hitTab) 

		missed.targets = readLines(paste('_missed.targets.', hiti, sep=''))
		n.missed.tar = length(missed.targets)
	}
	targets = missed.targets
}	

# final summary
writeLines(missed.targets, paste('missed.targets.cpu', sep=''))
system("cat _paired.ranked.* > paired.ranked.cpu")

system(paste('rm ../_running_', cpu, sep=''))

## THE END
