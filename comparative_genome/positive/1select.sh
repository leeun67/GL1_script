
for i in `cat 1`
do
	cat condem.null.txt |sed  "s/^      seqfile =/      seqfile = \.\/pp\/${i}.phy/g" |sed  "s/^      outfile =/      outfile = \.\/${i}.null.out/g" >${i}.null
	
	cat condem.alt.txt |sed  "s/^      seqfile =/      seqfile = \.\/pp\/${i}.phy/g" |sed  "s/^      outfile =/      outfile = \.\/${i}.alt.out/g" >${i}.alt

	codeml ${i}.null
	codeml	${i}.alt

done
