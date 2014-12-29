#
# This program calls primer3 to design primers (Right/reversed primer) using 
#   sequence templates in seq/ folder. This program is to be called by main.R.
# 
# Requires primer3path (should have been set by its master program)
# 
# Output: 
#   Primer candidates are output in out/ folder.


if (!exists('ncpu')){ncpu=8}
maxcpu=ncpu

system('rm -rf out')
system('mkdir out')

# cleanup
system('rm _running_* primer3.exome.setting* missed.seq_*')
nrun = 1

candidateDesign = function(design.round='initial', missedSeq=0
			   , seqFolder='seq'
			   , setting='primer3.exome.setting.1'){

	# get seqFile 
	if (design.round=='initial'){seqFile = dir(seqFolder)} ## initial design 
	if (design.round=='iteration'){seqFile = missedSeq} ## "missedSeq" is a resevered name
	if (design.round=='iteration' && missedSeq==0){stop} ##
	
	for (seq in seqFile){
		# max targets: maccpu
		cat('nrun: ', nrun, '\n')
		while (nrun  >= maxcpu){
			Sys.sleep(2)
			nrun=length(list.files('./',' _running_'))
		}

		runFile =  paste('_running_', seq, sep='')

		## jobs
		# prepare full paths for primer3 (I got errors
		#   when using reletive path)
		# if setting file not exists in current folder, make one
		setting.file = list.files(ampdir, setting, full=T)

		if (!file.exists(setting)){
			system(paste('grep -ve "^=" -ve PRIMER_MISPRIMING_LIBRARY -ve PRIMER_THERMODYNAMIC_PARAMETERS_PATH ', setting.file, " > _tmp.setting", sep=''))
			mispriming.file = list.files(depdir, NGSadaptors_and_humanRep, full=T)
			cat(paste("PRIMER_MISPRIMING_LIBRARY=", mispriming.file, sep='')
			    , file='_tmp.setting', append=T, sep='\n')
			primer3config = paste(list.dirs(primer3path, full=T)[1]
					      , '/src/primer3_config/', sep='')
			cat(paste("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=", primer3config
				  , sep=''), file='_tmp.setting', append=T, sep='\n')
			cat("=", file='_tmp.setting', append=T)
			system(paste('mv _tmp.setting ', setting, sep=''))
		}

		writeLines(paste(primer3path, "/src/primer3_core -p3_settings_file=", setting, " < ", seqFolder, "/", seq, " > out/", seq, "\n rm ", runFile, sep=''), runFile)
		system(paste('bash ', runFile, ' &'))

		nrun=length(list.files('./', '_running_'))
		cat(seq, '\n')
	}

		# wait
		while (nrun >0){
			Sys.sleep(5)
			nrun=length(list.files('./', '_running_'))
		}


	## if less than 10 candidates, send to missedSeq
	missedSeq = system("grep ok out/* | sed 's/:.* / /' | sed 's/out.//' | awk '{if ($2 < 10) print $1}'", intern=T)
	writeLines(missedSeq, paste('missed.seq', seqFolder, setting, sep='_'))

	(miss.n = length(missedSeq))
	(target.n = length(seqFile))
	(miss.rate = miss.n/target.n)
	}



####################
# call primer3 - round.1
candidateDesign(design.round='initial', missedSeq=0
		, seqFolder = 'seq'
		, setting='primer3.exome.setting.1')

	## if less than 10 candidates, send to missedSeq
	missedSeq = system("grep ok out/* | sed 's/:.* / /' | sed 's/out.//' | awk '{if ($2 < 10) print $1}'", intern=T)
	(miss.n = length(missedSeq))


# round 2
if (miss.n > 0){candidateDesign(design.round='iteration'
				, missedSeq = missedSeq
				, seqFolder='seq'
				, setting='primer3.exome.setting.2')}

	## if less than 10 candidates, send to missedSeq
	missedSeq = system("grep ok out/* | sed 's/:.* / /' | sed 's/out.//' | awk '{if ($2 < 10) print $1}'", intern=T)
	(miss.n = length(missedSeq))


# round 3
if (miss.n > 0){candidateDesign(design.round='iteration'
				, missedSeq = missedSeq
				, seqFolder='seq'
				, setting='primer3.exome.setting.3')}

	## if less than 10 candidates, send to missedSeq
	missedSeq = system("grep ok out/* | sed 's/:.* / /' | sed 's/out.//' | awk '{if ($2 < 10) print $1}'", intern=T)
	(miss.n = length(missedSeq))


# round 4 - noMask seq, setting 1
if (miss.n > 0){candidateDesign(design.round='iteration'
				, missedSeq = missedSeq
				, seqFolder='seq.noMask'
				, setting='primer3.exome.setting.1')}

	## if less than 10 candidates, send to missedSeq
	missedSeq = system("grep ok out/* | sed 's/:.* / /' | sed 's/out.//' | awk '{if ($2 < 10) print $1}'", intern=T)
	(miss.n = length(missedSeq))


# round 5 - noMask seq, setting 2
if (miss.n > 0){candidateDesign(design.round='iteration'
				, missedSeq = missedSeq
				, seqFolder='seq.noMask'
				, setting='primer3.exome.setting.2')}

	## if less than 10 candidates, send to missedSeq
	missedSeq = system("grep ok out/* | sed 's/:.* / /' | sed 's/out.//' | awk '{if ($2 < 10) print $1}'", intern=T)
	(miss.n = length(missedSeq))


# round 6 - noMask seq, setting 3
if (miss.n > 0){candidateDesign(design.round='iteration'
				, missedSeq = missedSeq
				, seqFolder='seq.noMask'
				, setting='primer3.exome.setting.3')}

	## if only 1 candidate, send to missedSeq
	missedSeq = system("grep -e 'ok 0' -e 'ok 1$' out/* | sed 's/:.*//' | sed 's/out\\///' ", intern=T)
	(miss.n = length(missedSeq))

## END of calling primer3
