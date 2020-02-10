#!/bin/bash 

#simple bash script for quantifying RNAseq data using Salmon.

curl -O https://www.aniseed.cnrs.fr/aniseed/download/?file=data%2Fcirobu%2FKHNCBI.Transcript.2018.fasta.zip
unzip \?file\=data%2Fcirobu%2FKHNCBI.Transcript.2018.fasta.zip
rm \?file\=data%2Fcirobu%2FKHNCBI.Transcript.2018.fasta.zip

salmon index -t KHNCBI.Transcript.2018.fasta -i Crob_index

for folder in AS05_Pax_rep*;
do
	for reads in $folder/AS05-*R1*.gz;
	do 
		salmon quant -i Crob_index -l A \
				 -1 ${reads} \
				 -2 ${reads/R1/R2} \
				 --gcBias \
				 --seqBias \
				 --numBootstraps 100 \
				 -p 8 --validateMappings -o /mnt/d/quants/${reads/_R1_001.fastq.gz/}_quant
	done
done 
