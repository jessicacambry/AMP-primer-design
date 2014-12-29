###########################
## pair GSP1-GPS2
##      by pair-panelty
###########################

currdir = getwd()
system('rm -rf pairing', ignore.stderr =T)
system('mkdir pairing')
setwd('pairing')

## generate list of targets from primer.psl
system("sort -k1,1 ../primer.psl > primer.psl.srt")
system("cut -d ' ' -f 1 primer.psl.srt > _targets")
system("sort -u _targets > target.list")
system("rm _targets")

## For multi-threading
pair.targets = readLines('target.list')
n.targets = length(pair.targets)
n.tar.cpu = ceiling(n.targets/ncpu)
seq1 = seq(1, n.targets, n.tar.cpu)
seq2 = seq(n.tar.cpu, n.targets, n.tar.cpu)
if (length(seq2) != length(seq1)){
    seq2 = c(seq2, n.targets)
}
n.seq = length(seq1)

# prep sub targets
for (cpu in 1:n.seq){
	system(paste('mkdir ', cpu, sep=''))
	target.n1 = seq1[cpu]
	target.n2 = seq2[cpu]
	targets.cpu = pair.targets[target.n1 : target.n2] 
	writeLines(targets.cpu, paste(cpu, '/targets', sep=''))
}

# multithread running
for (cpu in 1:n.seq){
        system(paste('touch _running_', cpu, sep=''))
	system(paste('Rscript --vanilla ', ampdir, '/scripts/pairing.GSPs.sub.R '
		     , cpu, ' ', ampdir, ' &', sep=''))
}
        #=== wait
	Sys.sleep(5)
        nrun=length(list.files('./', '_running_'))
        while (nrun >0){
                Sys.sleep(5)
                nrun=length(list.files('./', '_running_'))
        }

        ##=== all cpus paired
        system('cat */missed.targets.cpu > ../unpaired.targets')
        system('cat */paired.ranked.cpu > paired.ranked.allcpus')
        system('head -n1 paired.ranked.allcpus > ../ranked.all.candidate.pairs')
        system('grep -v penalty paired.ranked.allcpus >> ../ranked.all.candidate.pairs')

setwd(currdir)

# remove intermediate data
system('rm -rf pairing')

