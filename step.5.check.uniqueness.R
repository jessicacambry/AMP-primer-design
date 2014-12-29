##
## This program is to be called by main.R to
##   check uniqueness of GSP1s and GSP2s tail 12 bases
##
## input:
## 	- pairs1
## output:
##      - pairs.keep (unique)
##

chkt12 = function(d = pairs1){
	d$r1.n = nchar(d$r1.seq)
	d$r2.n = nchar(d$r2.seq)
	d$r1.t12 = substr(d$r1.seq, d$r1.n-11, d$r1.n)
	d$r2.t12 = substr(d$r2.seq, d$r2.n-11, d$r2.n)

	all.t12 = c(d$r1.t12, d$r2.t12)
	tab.t12 = data.frame('freq' = sort(table(all.t12)))
	tab.t12$t12 = row.names(tab.t12)
	tab.t12b = subset(tab.t12, freq >1)

	# not unique
	nu.t12s = tab.t12b$t12
	nutars = d$target[d$r1.t12 %in% nu.t12s | d$r2.t12 %in% nu.t12s]
	return(nutars)
}

(nutars = chkt12(d=pairs1))
n.new.nutar = length(nutars)

if (length(nutars) != 0){
	pairs.keep = subset(pairs1, !(target %in% nutars))

	# replace one target at a time, chk again
        for (nutar in nutars){
		while (n.new.nutar >0){
                system("head -n1 ranked.all.candidate.pairs > _tmp.cand ")
                system(paste("grep ", nutar, " ranked.all.candidate.pairs >> _tmp.cand"
                             , sep=''))
		tmpp = read.table('_tmp.cand', header=T, stringsAsFactors=F)
		tmpp$target = sub(":.*", "", tmpp$r1.qName)	

		# replace and chk again
		n.new.nutar = 1
		for (r in 2:nrow(tmpp)){
			if (n.new.nutar >0){
			print(r)
			Sys.sleep(1)
				pairs.new = rbind(pairs.keep, tmpp[r, ])
				new.nutars = chkt12(d=pairs.new)
				n.new.nutar = length(new.nutars)
				print(new.nutars)
				}
			}

		pairs.keep = pairs.new
			nutar.remain = nutars[-which(nutars == nutar)]
			pairs.remain = subset(pairs1, target %in% nutar.remain)
			pairs.keep.plus = rbind(pairs.keep, pairs.remain)
			new.nutars = chkt12(d=pairs.keep.plus)
			n.new.nutar = length(new.nutars)
			if (n.new.nutar == 0){
				pairs.keep = pairs.keep.plus
			}
		}
	} 
} else {
	pairs.keep = pairs1
}


## keep t12 for future
 pairs.keep$r1.n = nchar(pairs.keep$r1.seq)
 pairs.keep$r2.n = nchar(pairs.keep$r2.seq)
 pairs.keep$r1.t12 = substr(pairs.keep$r1.seq, pairs.keep$r1.n-11, pairs.keep$r1.n)
 pairs.keep$r2.t12 = substr(pairs.keep$r2.seq, pairs.keep$r2.n-11, pairs.keep$r2.n)

 all.t12 = c(pairs.keep$r1.t12, pairs.keep$r2.t12)

# End

