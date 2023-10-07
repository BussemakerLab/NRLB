#!/bin/bash
cd /Users/chaitanya/Documents/Research/SELEX/MultinomialPaper/ModelComparisons/dbsampling

startIdx=$(($1+1))
endIdx=$(($2+1))

baseDir=$(pwd)"/"

Emin=-15
Emax=15
divSize=.02
nInit=1000
fTol=1e-4
minCount=100
fCrit=.9

#Info file parsing
cat dbIDs.tsv | head -n+$endIdx | tail -n+$startIdx | while read p
do
        protName=$(echo $p | cut -f1 -d\ )
        round=$(echo $p | cut -f7 -d\ )
        mode=$(echo $p | cut -f8 -d\ )
        modelPath="../DeepBind/AUROC/FittedModels/"$protName"/R"$round"_M"$mode".dat"
        modelIndex=$(($(echo $p | cut -f9 -d\ )-1))
        nBases=$(echo $p | cut -f3 -d\ )
        varLen=$(echo $p | cut -f2 -d\ )
        lFlank=$(echo $p | cut -f4 -d\ )
        rFlank=$(echo $p | cut -f5 -d\ )
        echo "Curently processing "$protName
        java -Xmx2000m -cp ./bin/bioconductor/selex.jar:./bin/src:. utils.WangLandauDOS $baseDir $modelPath $modelIndex $nBases $varLen $lFlank $rFlank $Emin $Emax $divSize $nInit $fTol $minCount $fCrit > NRLBDensities/$protName.tsv
done
