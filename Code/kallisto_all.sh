#!/bin/bash

declare -a arr=("1" "2" "3" "4" "5" "6")
#declare -a arr=("7" "8" "9" "10" "11" "12")
#
#for i in "${arr[@]}"; do
#    echo ${i}
#    echo "Starting kallisto mapping to genome";
#    ~/workspace/kallisto/kallisto quant -i ~/workspace/kallisto/indices/neisseria_genome.idx -o ~/workspace/kallisto/neisseria_04_22_16/results/SQ${i} --single -l 50 -s 10 -b 2 --pseudobam ~/Projects/Neisseria/Data/SequencingData/Neisseria/FASTQ/SQ-${i}_no_rrnas.fastq > ~/workspace/kallisto/neisseria_04_22_16/SQ${i}.sam
#    
#    echo "Running my python script";
#    python sam_to_wig.py ~/workspace/kallisto/neisseria_04_22_16/SQ${i}.sam  
#    
#    echo "removing the SAM file";
#    rm ~/workspace/kallisto/neisseria_04_22_16/SQ${i}.sam
#    
#    echo "Starting kallisto mapping to CDS"
#    ~/workspace/kallisto/kallisto quant -i ~/workspace/kallisto/indices/neisseria_transcriptome.idx -o ~/workspace/kallisto/neisseria_04_22_16/results/SQ${i} --single -l 50 -s 10 -b 100 ~/Projects/Neisseria/Data/SequencingData/Neisseria/FASTQ/SQ-${i}_no_rrnas.fastq
#    
#done
#
#
#Testing potential transcriptome
#for i in "${arr[@]}"; do
#    echo ${i}
#    echo "Starting kallisto mapping to CDS"
#    ~/workspace/kallisto/kallisto quant -i ~/workspace/kallisto/indices/neisseria_potential_transcriptome.idx -o ~/workspace/kallisto/neisseria_05_02_16/results/SQ${i} --single -l 50 -s 10 -b 100 ~/Projects/Neisseria/Data/SequencingData/Neisseria/FASTQ/SQ-${i}_no_rrnas.fastq
#    
#done
##################################################
##################################################
#Re-doing normal transcriptome
#for i in "${arr[@]}"; do
#    echo ${i}
#    echo "Starting kallisto mapping to the genome"
#    ~/workspace/kallisto/kallisto quant -i ~/workspace/kallisto/indices/neisseria_genome.idx -o ~/workspace/kallisto/neisseria_02_28_17/results/SQ${i} --bias --single -l 50 -s 10 -b 2 -t 2 ~/Projects/Neisseria/Data/SequencingData/Neisseria/FASTQ/SQ-${i}_no_rrnas.fastq
#    echo "Starting kallisto mapping to CDS"
#    ~/workspace/kallisto/kallisto quant -i ~/workspace/kallisto/indices/neisseria_transcriptome.idx -o ~/workspace/kallisto/neisseria_02_28_17/results/SQ${i} --bias --single -l 50 -s 10 -b 2 -t 2 ~/Projects/Neisseria/Data/SequencingData/Neisseria/FASTQ/SQ-${i}_no_rrnas.fastq
#
#done
##################################################
##################################################
#SAM files of transcriptome output for visualization of paralogs
for i in "${arr[@]}"; do
    echo ${i}
    echo "Starting kallisto mapping to CDS"
    ~/workspace/kallisto/kallisto quant -i ~/workspace/kallisto/indices/neisseria_transcriptome.idx -o ~/workspace/kallisto/neisseria_04_20_17/results/SQ${i} --bias --single -l 50 -s 10 -b 2 -t 2 --pseudobam ~/Projects/Neisseria/Data/SequencingData/Neisseria/FASTQ/SQ-${i}_no_rrnas.fastq > ~/workspace/kallisto/neisseria_04_20_17/SQ${i}.sam

done
