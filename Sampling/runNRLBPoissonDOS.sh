#!/bin/bash

startIdx=$(($1+1))
endIdx=$(($2+1))

baseDir="/Users/chaitanya/Documents/Research/SELEX/MultinomialPaper/ModelComparisons/"

Emin=-30
Emax=30
divSize=.02
nInit=1000
fTol=1e-4
minCount=100
fCrit=.9

#Info file parsing
cat /Users/chaitanya/Desktop/dbsampling/dbIDs.tsv | head -n+$endIdx | tail -n+$startIdx | while read p
do
	protName=$(echo $p | cut -f1 -d\ )
	round=$(echo $p | cut -f7 -d\ )
	mode=$(echo $p | cut -f8 -d\ )
	modelPath="DeepBind/AUROC/FittedModels/"$protName"/R"$round"_M"$mode".dat"
	R0FName=$(ls "/Users/chaitanya/Documents/Research/SELEX/MultinomialPaper/ModelComparisons/DeepBind/AUROC/FittedModels/"$protName | grep "R0.dat")
	R0Path="DeepBind/AUROC/FittedModels/"$protName"/"$R0FName
	modelIndex=$(($(echo $p | cut -f9 -d\ )-1))
	probesPath="../../../../../Desktop/dbsampling/CountScores/"$protName"_R0R1.tsv"
	echo "Curently processing "$protName
	java -Xmx2000m -cp ./bin/bioconductor/selex.jar:./bin/src:. utils.WangLandauDOSProbeScore $baseDir $modelPath $modelIndex $R0Path $Emin $Emax $divSize $nInit $fTol $minCount $fCrit $probesPath > NRLBPoissonDensities/$protName.tsv
done
