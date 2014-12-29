
cat split/psl.* > __all.psl

## psl file format:
# $1 matches - Number of matching bases that aren't repeats.
# $11 qSize - Query sequence size.
# $13 qEnd - Alignment end position in query.
# So, to keep only alignments that matters - with complete match 
#   of 12 bases at the 3' end:

awk '{if ($1>=12 && $11==$13) print}' __all.psl > __all.psl.matt
sort -k10,10b -k1,1nr -k2,2n -k20,20n __all.psl.matt \
	| sed 's:_SEQUENCE::' | sort -k10,10b > all.psl.matt.sorted
rm __*

grep EXPLAIN  out/* |  wc              #
grep EXPLAIN  out/* | grep 'ok 0' | wc #
grep EXPLAIN  out/* | grep -v 'ok 0' | wc #  see if all primers were hit

grep SEQUENCE= out/* -A 1 | grep ',' > __primer.pos
sed 's/,/ /' __primer.pos | sed 's:out/::' | sed 's/=/ /' \
	| sed 's/-PRIMER/:PRIMER/' > __primer.pos2

# TM
grep _TM= out/* > __primer.tm
sed 's/,/ /' < __primer.tm | sed 's:out/::' | sed 's/_TM=/ /' > __primer.tm2

# mispriming
grep _MISPRIMING= out/* | tr "'" "-" > __primer.mis
sed 's/,/ /' __primer.mis | sed 's:out/::' | sed 's/  / /' \
	| sed 's/  / /' | sed 's/ /zzz/' | sed 's/_LIBRARY_.*zzz/zzz/' \
	| sed 's/reverse //' | sed 's/ .*//' | sed 's/zzz/ /' > __primer.mis2

# seq
grep _SEQUENCE= out/* > __primer.seq
sed 's/,/ /' __primer.seq | sed 's:out/::' | sed 's/_SEQUENCE=/ /' > __primer.seq2

sort -k1b,1 __primer.pos2 > __primer.pos3
sort -k1b,1 __primer.tm2 > __primer.tm3
sort -k1b,1 __primer.mis2 > __primer.mis3
sort -k1b,1 __primer.seq2 > __primer.seq3

join __primer.pos3 __primer.tm3 > __primer_1
join __primer_1 __primer.mis3 > __primer_2
join __primer_2 __primer.seq3 > __primer_3

sort -k1b,1 __primer_3 > primer.pos.tm.sorted
join -2 10 primer.pos.tm.sorted all.psl.matt.sorted > __psl.primer3
sed 's/:/ /' __psl.primer3 | awk '{posLen=$3","$4; print $0,posLen}' > primer.psl

rm __*

