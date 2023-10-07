#!/bin/bash
#for the cluster

tail -n+2 YangDataInfo.txt | while read p
do
	protName=$(echo $p | cut -f1 -d\ )
	R0File=$(echo $p | cut -f3 -d\ )
	R1File=$(echo $p | cut -f4 -d\ )
	if [[ $R1File = NA ]]
        then
                echo "R1 file for "$protName" is missing. Skipping..."
        else
		echo "Copying files for "$protName"..."
                cp /vega/hblab/users/htr2104/multiRoundModel/rawData/Yang2017/$R0File .
		cp /vega/hblab/users/htr2104/multiRoundModel/rawData/Yang2017/$R1File .
		echo "Finished copying files. Computing count tables..."
		./buildTable.py --gz --fastq $R0File $R1File > $protName"_R0R1.tsv"
		echo "Counting complete."
		rm $R0File $R1File
        fi
done
ls *_R0R1* | while read p; do sed -i '/N/d' $p; done;
echo "Data pull complete."
