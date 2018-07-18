#!/bin/bash
#Ensure the proper number of config files
if [[ $# -lt 2 ]] || [[ $# -gt 3 ]]
then 
	echo "Please provide both data and NRLB configuration files!"
	exit 1
fi

if [[ $# -eq 3 ]]
then
	memUsage="Xmx$3"
else
	memUsage="Xmx12000"
fi
memUsage+="m"

#Ensure OS Compatability
replace() {
	if [[ $(uname) != Linux ]]
	then
     		sed -i '' "$@"
	else
        	sed -i -e "$@"
	fi
}

homeDir=$(pwd)
dataConfig=$1
nrlbConfig=$2
protName=$(grep "protName" $dataConfig | cut -f2 -d=)
tmpXMLFileName="tempXML_$protName.xml"
tmpConfigName="tempNRLB_$protName.config"
forceRecalc=$(grep "forceRecalc" $dataConfig | cut -f2 -d=)
#clear
cd $homeDir

#Copy the temporary xml file in the tmp folder
mkdir tmp 2> /dev/null
cp bin/template.xml tmp/$tmpXMLFileName
#Replace necessary elements in xml file
grep "=" $dataConfig  | while read line
do 
	var1=$(echo $line | cut -f1 -d= | xargs)
	var2=$(echo $line | cut -f2 -d= | xargs)
	replace "s~$var1~$var2~g" tmp/$tmpXMLFileName
done
#Run kmer counting for both R0 and R1
echo "Counting Kmers..."
varLen=$(grep "variableRegionLength" $dataConfig | cut -f2 -d=)
#Count kmers only if the file does not exist. Count R0 File first
if [ ! -e tmp/$protName.0.R0.0.$varLen.dat ] || [[ $forceRecalc = true ]]
then
	java -Xmx8000m -cp ./bin/NRLB.jar utils.Build_MM_Kmer_Table $varLen 0 $homeDir/tmp/ $homeDir/tmp/$tmpXMLFileName $protName.0 R0 > /dev/null
#	java -Xmx8000m -cp ./bin/bioconductor/selex.jar:./bin/src/:. utils.Build_MM_Kmer_Table $varLen 0 $homeDir/tmp/ $homeDir/tmp/$tmpXMLFileName $protName.0 R0 > /dev/null
fi
#STATUS="${?}"
if [ $? -ne 0 ]
then
	exit $?
fi
#Count Kmers for R1 File
if [ ! -e tmp/$protName.1.R1.1.$varLen.dat ] || [[ $forceRecalc = true ]]
then
	java -Xmx8000m -cp ./bin/NRLB.jar utils.Build_MM_Kmer_Table $varLen 1 $homeDir/tmp/ $homeDir/tmp/$tmpXMLFileName $protName.1 R1 > /dev/null
#	java -Xmx8000m -cp ./bin/bioconductor/selex.jar:./bin/src/:. utils.Build_MM_Kmer_Table $varLen 1 $homeDir/tmp/ $homeDir/tmp/$tmpXMLFileName $protName.1 R1 > /dev/null 
fi
#STATUS="${?}"
if [ $? -ne 0 ]
then
	exit $?
fi
#Remove temp xml file and other temporary files
rm tmp/*.xml 2> /dev/null
rm -rf tmp/log 2> /dev/null
rm tmp/*.gz-* 2> /dev/null
echo "Counting complete. Constructing Round0 Model..."

#Begin constructing R0 Model
mkdir R0Models 2> /dev/null
lFlank=$(grep "lFlank" $dataConfig | cut -f2 -d=)
rFlank=$(grep "rFlank" $dataConfig | cut -f2 -d=)
R0MinK=$(grep "R0MinK" $dataConfig | cut -f2 -d= | tr -d '\040\011\012\015')
R0MaxK=$(grep "R0MaxK" $dataConfig | cut -f2 -d= | tr -d '\040\011\012\015')
if [ $R0MaxK -gt 0 ] 2> /dev/null
then 
	R0MK=$R0MaxK
else 
	R0MK="null"
fi
if [ ! -e R0Models/$protName.dat ] || [[ $forceRecalc = true ]]
then
	rm R0Models/$protName.* 2> /dev/null
	if [[ $R0MK = null ]]
	then
		java -Xmx8000m -cp ./bin/NRLB.jar model.Round0Regression $varLen $homeDir/tmp/ $protName $lFlank $rFlank $homeDir/R0Models/ $R0MinK > /dev/null
	else
		java -Xmx8000m -cp ./bin/NRLB.jar model.Round0Regression $varLen $homeDir/tmp/ $protName $lFlank $rFlank $homeDir/R0Models/ $R0MinK $R0MK > /dev/null
	fi
fi
if [ $? -ne 0 ]
then 
	exit $?
fi
echo "Round0 Model Construction Complete. Starting NRLB..."

#Begin NRLB Model Construction
forceRecalc=$(grep "forceRecalc" $nrlbConfig | cut -f2 -d=)
#Copy temporary configuration file and fill in relevant information
mkdir NRLBModels 2> /dev/null
cp bin/template.config NRLBModels/$tmpConfigName

replace "s~@workingDir~$homeDir/tmp/~g" NRLBModels/$tmpConfigName
replace "s~@seqRunName~$protName.1~g" NRLBModels/$tmpConfigName
replace "s~@sampleName~R1.1~g" NRLBModels/$tmpConfigName
replace "s~@varLen~$varLen~g" NRLBModels/$tmpConfigName
replace "s~@lFlank~$lFlank~g" NRLBModels/$tmpConfigName
replace "s~@rFlank~$rFlank~g" NRLBModels/$tmpConfigName
replace "s~@R0ModelPath~$homeDir/R0Models/$protName.dat~g" NRLBModels/$tmpConfigName
replace "s~@outputPath~$homeDir/NRLBModels/~g" NRLBModels/$tmpConfigName
toutName=${nrlbConfig##*/}
nrlbOutName=${toutName%.*}
replace "s~@outputName~$protName-$nrlbOutName~g" NRLBModels/$tmpConfigName

#Replace necessary elements in config file
grep "=" $nrlbConfig  | while read line
do
        var1=$(echo $line | cut -f1 -d= | xargs)
        var2=$(echo $line | cut -f2 -d= | xargs)
        replace "s~@$var1~$var2~g" NRLBModels/$tmpConfigName
done
cat NRLBModels/$tmpConfigName | while IFS= read i; do if [[ $i = @MSLOC ]]; then cat $nrlbConfig | grep "modeSeed"; else echo "$i"; fi; done > tmp.txt
mv tmp.txt NRLBModels/$tmpConfigName

#Run NRLB
if [ ! -e NRLBModels/$protName-$nrlbOutName.csv ] || [[ $forceRecalc = true ]]
then
	rm NRLBModels/$protName-$nrlbOutName.* 2> /dev/null
	java -$memUsage -cp ./bin/NRLB.jar model.NRLB NRLBModels/$tmpConfigName
fi

#Remove config file
rm NRLBModels/$tmpConfigName
