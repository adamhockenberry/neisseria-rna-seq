from Bio import SeqIO
import argparse
import numpy as np

if __name__ == '__main__':
    ######Parse input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    args = parser.parse_args()
    # print(args.filename)
    assert args.filename.count('.sam') == 1
    #######Set free parameters
    length_restriction = (20, 100)
    #######Get file names to write
    wiggle_out_f = args.filename.replace('.sam', '_f.wig')
    wiggle_out_r = args.filename.replace('.sam', '_r.wig')
    wiggle_out_normed_f = args.filename.replace('.sam', 'normed_f.wig')
    wiggle_out_normed_r = args.filename.replace('.sam', 'normed_r.wig')
    #######Read the genome
    genome_list = list(SeqIO.parse('../Data_release/fa1090.gb', 'genbank'))
    assert len(genome_list) == 1
    genome = genome_list[0]
    #######Instantiate empty dictionaries

    genome_length = len(genome.seq)
    fwd_array = np.zeros(genome_length, dtype=np.int) #These arrays will contain position:read_count pairs for all mapped locations
    rev_array = np.zeros(genome_length, dtype=np.int)

    with open(args.filename, 'r') as infile:
        for line in infile:
            if line[0] == '@': #Ignore these lines
                continue
            split_line = line.split('\t')
            if split_line[1] == '4': #Ignore unmapped reads
                continue
            seq = split_line[9] #In the sam file position 10 (index 9) should contain the mapped sequence
            seq_length = len(seq)
            if seq_length < length_restriction[0] or seq_length > length_restriction[1]: #ignore reads outside of length restriction
                continue
            #####Get positions of the reads
            start_loc = int(split_line[3])-1 #I noticed the subtraction of one is necessary for proper mapping
            end_loc = start_loc+seq_length
            #####Add reads to relevant dictionaries
            if split_line[1] == '0':
                fwd_array[start_loc:end_loc] += 1
                if end_loc >= genome_length:
                    fwd_array[0:end_loc%genome_length] += 1
            elif split_line[1] == '16':
                rev_array[start_loc:end_loc] += 1
                if end_loc >= genome_length:
                    rev_array[0:end_loc%genome_length] += 1

    total_mapped = 0
    for feature in genome.features:
        if feature.type == 'CDS':
            start, end = feature.location.start, feature.location.end
            if feature.strand == 1:
                total_mapped += sum(fwd_array[start:end])
            elif feature.strand == -1: 
                total_mapped += sum(rev_array[start:end])
    
    norm_constant = total_mapped / 1000000000.

    fwd_array_normed = fwd_array / norm_constant
    rev_array_normed = rev_array / norm_constant
    
    with open(wiggle_out_f, 'w') as f, \
            open(wiggle_out_normed_f, 'w') as fn, \
            open(wiggle_out_r, 'w') as r, \
            open(wiggle_out_normed_r, 'w') as rn:
        f.write('variableStep chrom=NC_002946 \n')
        fn.write('variableStep chrom=NC_002946 \n')
        r.write('variableStep chrom=NC_002946 \n')
        rn.write('variableStep chrom=NC_002946 \n')
        for i, (fv, fvn, rv, rvn) in enumerate(zip(fwd_array, fwd_array_normed, rev_array, rev_array_normed)):
            if fv:
                f.write('{}\t{}\n'.format(i, fv))
                fn.write('{}\t{}\n'.format(i, fvn))
            if rv:
                r.write('{}\t{}\n'.format(i, rv))
                rn.write('{}\t{}\n'.format(i, rvn))
