
## Map candidate primers against reference genome using BLAT.

## For BLAT manual, check out https://users.soe.ucsc.edu/~kent/src/

## parameters set by master program
# t = 11
# s = 4
# n.srv = ceiling(ncpu/4)
# ClnPsrv = 4
# n.fa=n.srv*ClnPsrv
# repMatch=1024
# start.port=9000
# faDir = 'split'
# depdir = '~/dependency-AMP'
# pjdir = '~/project-AMP'
# db = 'hg19.2bit'

## Initiate gfServer if not already up
currdir = getwd()
setwd(pjdir)

allports = seq(start.port, start.port + n.srv -1)
n.ports.to.setup = n.srv
ports.being.setup = NA

while (n.ports.to.setup !=0){
	ports.to.setup = allports

	for (port in ports.to.setup){
		# check status
		port.status = system(paste(blatdir, "/gfServer status localhost "
				     , port, sep=''), intern=T, ignore.stderr=T)
		# setup if not ready AND not being setup
		if (length(port.status) == 0 
		     & !(port %in% ports.being.setup)
		    ){
                      system(paste(blatdir, "/gfServer start localhost "
				   , port, " ", depdir, "/", db, " -tileSize=", t
				   , " -stepSize=", s, " -repMatch=", repMatch
				   , " -canStop &", sep=''))
			ports.being.setup = c(ports.being.setup, port)
		      } else {
		           ports.to.setup = ports.to.setup[-which(ports.to.setup==port)]
		      }
	}

	n.ports.to.setup = length(ports.to.setup) 
	if (n.ports.to.setup != 0){
		Sys.sleep(120)
	}
}

ports.being.setup = ports.being.setup[-1] ## remove NA

## wait until all set
allset = 0
while (allset == 0){

	allset =1
	# if any of ports not done, let allset = 0
        for (port in ports.being.setup){
	       port.status = system(paste(blatdir, "/gfServer status localhost "		                                        , port, sep=''), intern=T, ignore.stderr=T)
	       if (length(port.status) == 0){
		       allset = 0
	       }
	}

	Sys.sleep(5)
}

setwd(currdir)


##
## split fasta
##

system('rm -rf split')
system('mkdir split')
system("grep SEQUENCE= out/* | sed 's/out\\//>/' | tr '=' '\n' > primer.candidates.0")
system(paste('R CMD BATCH --no-restore --no-save ', '"', "--args n.fa=", n.fa," inFa='primer.candidates.0' outDir='", faDir, "'", '" ', ampdir, '/scripts/splitFa.R', sep=''))

## BLAT
system('rm _job_*', ignore.stderr=TRUE)

# minIden = 93
# minScore = 12
# NAtype = 'dna'
# maxIntron = 700001

for (i in 0:(n.fa-1)){
        ij = 1000 + i
        port = 9000 + floor(i/ClnPsrv)
	runFile = paste('_job_', i, sep='_')

	writeLines(paste(blatdir, "/gfClient localhost ", port, " '/' ", faDir, "/fa."
			 , ij, " -t=", NAtype, " -out=psl -nohead -minScore=", minScore
			 , " -minIdentity=", minIden, " -maxIntron=", maxIntron, " "
			 , faDir, "/psl.", ij, " \n rm ", runFile, sep=''), runFile)
        system(paste('bash ', runFile, ' &'))
}

    	#======== wait
	nrun=length(list.files('./', '_job_'))
                while (nrun >0){
                        Sys.sleep(5)
                        nrun=length(list.files('./', '_job_'))
                }

if (keep_gfSvr != 1){
	system('killall -9 gfServer')
}

