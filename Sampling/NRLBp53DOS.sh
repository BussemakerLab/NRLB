#!/bin/bash

Emin=-13
Emax=13
divSize=.005
nInit=1000
fTol=1e-4
minCount=100
fCrit=.9
        
baseDir="/Users/chaitanya/Desktop/"
modelPath="p53-unmeth-WT.dat"
modelIndex=25
nBases=-62
varLen=26
lFlank="GT"
rFlank="AC"

java -Xmx4000m -cp ./bin/bioconductor/selex.jar:./bin/src:. utils.WangLandauDOS $baseDir $modelPath $modelIndex $nBases $varLen $lFlank $rFlank $Emin $Emax $divSize $nInit $fTol $minCount $fCrit > $baseDir/p53dist.tsv
