#!/bin/bash

cd /Users/chaitanya/Documents/Research/SELEX/MultinomialPaper/ModelComparisons/dbsampling

startIdx=$(($1+1))
endIdx=$(($2+1))

Emin=-40
Emax=40
divSize=1
nInit=1000
fTol=1e-4
minCount=100
fCrit=.9

#Info file parsing
cat dbIDs.tsv | head -n+$endIdx | tail -n+$startIdx | while read p
do
	protName=$(echo $p | cut -f1 -d\ )
	dbID=$(echo $p | cut -f6 -d\ )
	varLen=$(echo $p | cut -f2 -d\ )
	lFlank=$(echo $p | cut -f4 -d\ )
	rFlank=$(echo $p | cut -f5 -d\ )
	echo "Curently processing "$protName
	./deepbind $dbID $varLen $lFlank $rFlank $Emin $Emax $divSize $nInit $fTol $minCount $fCrit > DBDensities/$protName.tsv
done
