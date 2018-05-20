from Bio import SeqIO
import argparse

def get_intragenic_seqs(genome, buffer_length=30):
    starting_indent = 0 
    
    intragenic_seqs = []
    all_seq_lengths = []
    
    half_genome_length = len(genome.seq)/2.
    segments = []
    for feature in genome.features:
       try:
           if feature.type == 'gene':
               start = feature.location.start
               stop = feature.location.end
               if stop-start > half_genome_length:
                   continue
               segments.append((start-buffer_length, stop+buffer_length))
               if stop < start:
                   starting_indent += 1
       except:
           pass

    indices = []
    for start, finish in segments:
        indices.append(start)
        indices.append(-finish)
    indices = sorted(indices, key=lambda x: abs(x))
    indices.append(len(str(genome.seq))-1)

    temp_seqs = []
    temp_seq_names = []
    previous = 0 
    indentation = starting_indent
    for index in indices:
        if indentation == 0:
            if index != previous:
                temp_seqs.append(str(genome.seq[previous:index]))
                temp_seq_names.append('{}...{}'.format(previous, index))
        if index >= 0:
            previous = index
            indentation += 1
        else:
            previous = -index
            indentation -= 1
    intragenic_seqs_dict = dict(zip(temp_seq_names, temp_seqs))
    return intragenic_seqs_dict

if __name__ == '__main__':
    ######Parse input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    args = parser.parse_args()
    print(args.filename, args.filename.count('.gb')==1)
    assert args.filename.count('.gb') == 1, print("Are you sure this file ends in .gb or .gbk?")
    file_ending = args.filename.split('.')[-1]
    
    
    #######Get file names to write
    genome_fasta_out = args.filename.replace('.{}'.format(file_ending), '_genome.fasta')
    transcriptome_fasta_out = args.filename.replace('.{}'.format(file_ending), '_transcripts.fasta')
    anti_transcriptome_fasta_out = args.filename.replace('.{}'.format(file_ending), '_anti_transcripts.fasta')
    possible_transcriptome_fasta_out = args.filename.replace('.{}'.format(file_ending), '_potential_transcripts.fasta')
    
    #######Read the genome
    genome_list = list(SeqIO.parse(args.filename, 'genbank'))
    assert len(genome_list) == 1, print("Code can currently only process .gb files with one sequence record")
    genome = genome_list[0]

    output_handle = open(genome_fasta_out, "w")
    output_handle.write(">" + genome.id + '\n')
    output_handle.write(str(genome.seq))
    output_handle.close()

    output_handle = open(transcriptome_fasta_out, "w")
    for feature in genome.features:
        if feature.type == 'CDS':
            if 'locus_tag' in feature.qualifiers:
                gene_name = feature.qualifiers['locus_tag'][0]
            else:
                continue
        else:
            continue
        output_handle.write(">" + gene_name + '\n')
        output_handle.write(str(feature.extract(genome.seq))+ '\n')
    output_handle.close()
    
    output_handle = open(anti_transcriptome_fasta_out, "w")
    for feature in genome.features:
        if feature.type == 'CDS':
            if 'locus_tag' in feature.qualifiers:
                gene_name = feature.qualifiers['locus_tag'][0]
            else:
                continue
        else:
            continue
        output_handle.write(">" + gene_name + '\n')
        output_handle.write(str(feature.extract(genome.seq).reverse_complement())+ '\n')
    output_handle.close()
    
    output_handle = open(possible_transcriptome_fasta_out, "w")
    intragenic_seqs_dict = get_intragenic_seqs(genome)
    for key, value in intragenic_seqs_dict.items():
        output_handle.write(">{}\n".format(key))
        output_handle.write('{}\n'.format(value))
    for feature in genome.features:
        if feature.type == 'CDS':
            if 'locus_tag' in feature.qualifiers:
                gene_name = feature.qualifiers['locus_tag'][0]
            else:
                continue
        else:
            continue
        output_handle.write(">" + gene_name + '\n')
        output_handle.write(str(feature.extract(genome.seq))+ '\n')
    output_handle.close()
        
