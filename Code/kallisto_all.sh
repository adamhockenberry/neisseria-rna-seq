#!/bin/bash


########
########
#Should be pretty obvious but this is the code I used to run kallisto. All directories will need to be re-named accordingly
#if you plan on running this in any meaningful way
#Keeping this file here, however, for posterity's sake lest anyone be curious
########
########

declare -a arr=("1" "2" "3" "4" "5" "6")
#declare -a arr=("7" "8" "9" "10" "11" "12")
#
for i in "${arr[@]}"; do
    echo ${i}
    echo "Starting kallisto mapping to genome";
    ~/workspace/kallisto/kallisto quant -i ~/workspace/kallisto/indices/neisseria_genome.idx -o ~/workspace/kallisto/neisseria_04_22_16/results/SQ${i} --single -l 50 -s 10 -b 2 --pseudobam ~/Projects/Neisseria/Data/SequencingData/Neisseria/FASTQ/SQ-${i}_no_rrnas.fastq > ~/workspace/kallisto/neisseria_04_22_16/SQ${i}.sam
    
    echo "Running my python script";
    python sam_to_wig.py ~/workspace/kallisto/neisseria_04_22_16/SQ${i}.sam  
    
    echo "removing the SAM file";
    rm ~/workspace/kallisto/neisseria_04_22_16/SQ${i}.sam
    
    echo "Starting kallisto mapping to CDS"
    ~/workspace/kallisto/kallisto quant -i ~/workspace/kallisto/indices/neisseria_transcriptome.idx -o ~/workspace/kallisto/neisseria_04_22_16/results/SQ${i} --single -l 50 -s 10 -b 100 ~/Projects/Neisseria/Data/SequencingData/Neisseria/FASTQ/SQ-${i}_no_rrnas.fastq
done
