##
## Pairing GSP1-GSP2 based on a pair penalty, for one target. Pairs
##   are first created using a fast matrix method.
##
## input:
##	good.hit.tar.psl - which is BLAT result with primers
##           that show multiple hits in genome filtered.
##      good.hit.tar     - name of the good target e.g. 'ALK_001_1621_AntiSense'
## output:
##	paired.cands.ranked_$target	

good.hit.tar.psl = dget('__good.hit.tar.psl')
good.hit.tar = good.hit.tar.psl$tar[1]

miss = good.hit.tar.psl$MIS # misspriming source (repetitive, NGSadaptors, etc)
len = length(miss)
miss.m = matrix(miss!=rep(miss, each=len), nrow=len)
ll = 1:(len*len)
ll_missF = ll[miss.m[TRUE]]

if (length(ll_missF) >= 1){

x = ceiling(ll_missF / len)
y = as.numeric(ll_missF %% len)
y[y==0] = len
heads = good.hit.tar.psl$head
pairData = data.frame('x' = x, 'y' = y)
pairData$pos.x = heads[x]
pairData$pos.y = heads[y]
pairData$pos.d = pairData$pos.x - pairData$pos.y

## requires at least 5 base distance between GSP1 and GPS2 3' ends 
##   for increased enrichment specificity
pairD2 = subset(pairData, pos.d >= 5)

if (nrow(pairD2) >= 1){
	pairD2$pos.d30 = abs(pairD2$pos.d - 30)
	pairD2 = pairD2[order(pairD2$pos.d30), ]
	pairD2$pos.dd = floor(pairD2$pos.d30 / 10)
	pairD2$pos.yy = floor(pairD2$pos.y / 10)
	pairD2$penalty = pairD2$pos.yy + pairD2$pos.dd

	pairD3 = pairD2[order(pairD2$penalty, -pairD2$y),]

	pairD3$r1.qName = good.hit.tar.psl$qName[pairD3$x]
	pairD3$r1.seq = good.hit.tar.psl$candSeq[pairD3$x]
	pairD3$r1.TM = good.hit.tar.psl$TM[pairD3$x]
	pairD3$r1.chr = good.hit.tar.psl$tName[pairD3$x]
	pairD3$r1.tStart = good.hit.tar.psl$tStart[pairD3$x]
	pairD3$r1.tEnd = good.hit.tar.psl$tEnd[pairD3$x]

	pairD3$r2.qName = good.hit.tar.psl$qName[pairD3$y]
	pairD3$r2.seq = good.hit.tar.psl$candSeq[pairD3$y]
	pairD3$r2.TM = good.hit.tar.psl$TM[pairD3$y]
	pairD3$r2.chr = good.hit.tar.psl$tName[pairD3$y]
	pairD3$r2.tStart = good.hit.tar.psl$tStart[pairD3$y]
	pairD3$r2.tEnd = good.hit.tar.psl$tEnd[pairD3$y]

	pairD3$TM.d65 = abs(round((67 - pairD3$r2.TM) / 2))
	pairD3$TM.d12 = (round(pairD3$r2.TM) - round(pairD3$r1.TM)) / 2 - 1
	pairD3$TM.d12p = ifelse(pairD3$TM.d12 < 0, 100, pairD3$TM.d12)

	  ### pair penalty
	paired.cands = pairD3[order(pairD3$penalty, pairD3$TM.d65, pairD3$TM.d12p
				    , -pairD3$TM.d12),]
	paired.cands$rank = order(1:nrow(paired.cands))

	write.table(paired.cands, paste('__paired.cands.ranked_', good.hit.tar, sep='')
		    , row.names=F, col.names=T, quote=F)

	rm(good.hit.tar)
	rm(good.hit.tar.psl)
}
}
