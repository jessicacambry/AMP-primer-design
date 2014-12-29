## split fasta file into small files of equal size
## requires:
##	n.fa
##	inFa
##	ourDir
##      ampdir

args=(commandArgs(TRUE))
for(i in 1:length(args)){eval(parse(text=(args[[i]])))}

## input file is:
cat('Input file: ', inFa,'\n')

n.reads = as.numeric(system(paste("wc -l ", inFa, " | cut -d ' ' -f 1", sep=''), intern=T))/2
cat("### Number of reads:\t\t n=", n.reads, "\n")
n.split =2*(ceiling(n.reads/n.fa))
system(paste('rm -rf ', outDir, sep=''))
system(paste('mkdir ', outDir, sep=''))
system(paste("split ", inFa," ", outDir, "/fa.1 -a 3 -d -l ", n.split, sep=''))
