#!/bin/bash

startIdx=$(($1+1))
endIdx=$(($2+1))

baseDir="/Users/chaitanya/Documents/Research/SELEX/MultinomialPaper/ModelComparisons/DeepBind/"
probesPath="../../../../../../Desktop/dbsampling/probes.txt"

# Loop over all proteins
head -n+$endIdx YangDataInfo.tsv | tail -n+$startIdx | while read p
do
	protName=$(echo $p | cut -f1 -d\ )
	echo $protName
	R1FName=$(echo $p | cut -f4 -d\ )
	if [[ $R1FName = NA ]]
	then
		echo "Skipping scoring of "$protName"; no R1 file exists."
	else 
		echo "Scoring "$protName"..."
		# Extract information
		lFlank=$(grep $protName dbIDs.tsv | cut -f4)
		rFlank=$(grep $protName dbIDs.tsv | cut -f5)
		dbID=$(grep $protName dbIDs.tsv | cut -f6)
		round=$(grep $protName dbIDs.tsv | cut -f7)
        	mode=$(grep $protName dbIDs.tsv | cut -f8)
        	modelPath="../DeepBind/AUROC/FittedModels/"$protName"/R"$round"_M"$mode".dat"
        	modelIndex=$(($(grep $protName dbIDs.tsv | cut -f9)-1))
    		nBases=$(grep $protName dbIDs.tsv | cut -f3)
		# Extract probes
		cat "R2ROC/"$protName"_ROCR2.tsv" | cut -f1 | sed "s/^/$lFlank/" | sed "s/$/$rFlank/" > probes.txt
		# Compute DB Scores
		/Users/chaitanya/Documents/Research/SELEX/MultinomialPaper/ModelComparisons/DeepBind/deepbind/./deepbind $dbID < probes.txt > DBscores.txt
		# Compute NRLB Scores
		java -Xmx2000m -cp ./bin/bioconductor/selex.jar:./bin/src:. utils.ModelProfiler2 $baseDir $modelPath $modelIndex $probesPath $nBases > NRLBscores.txt
		# Merge
		paste DBscores.txt NRLBscores.txt | sed "1s/.*/DeepBind$(printf '\t')NRLB/" > "R2ROC/"$protName"_ROC_scores.tsv"
		rm DBscores.txt NRLBscores.txt probes.txt
	fi
done
echo "Finished scoring all probes."
