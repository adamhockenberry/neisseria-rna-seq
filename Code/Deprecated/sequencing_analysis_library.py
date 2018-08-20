import numpy as np
from matplotlib import pyplot as plt
import datetime
import os
from scipy import stats

plt.rcParams['axes.color_cycle'] = ['b', 'r', 'c', 'y', 'm', 'g', 'k']


def sam_to_wiggle(sam_file_name, organism, genome_list, length_restriction=(15,100), assignment='all'):
    '''
    This function reads in a .sam file, extracts the sequencing information and writes it to two
    separate .wig files (fwd and rev strands) as well as a summary.txt file
    '''
    fwd_dicty = {} #These dictionaries will contain position:read_count pairs for all mapped locations
    rev_dicty = {}
    mapped_len_dict = {}
    mapped_len_dict['Plus'] = []
    mapped_len_dict['Minus'] = []
    mapped_len_dict['Unmapped'] = []
    infile = open(sam_file_name)
    for line in infile:
        if line[0] == '@':
            continue
        split_line = line.split('\t')
        seq = split_line[9] #In the sam file position 10 (index 9) should contain the mapped sequence
        seq_length = len(seq)
        if split_line[1] == '4':
            mapped_len_dict['Unmapped'].append(seq_length)
            continue
        if seq_length < length_restriction[0] or seq_length > length_restriction[1]:
            mapped_len_dict['Unmapped'].append(seq_length)
            continue
        start_loc = int(split_line[3])-1 #I noticed the subtraction of one is necessary for proper mapping
        end_loc = start_loc+seq_length
        contributions = 1./seq_length
        if split_line[1] == '0':
            mapped_len_dict['Plus'].append(seq_length)
        elif split_line[1] == '16':
            mapped_len_dict['Minus'].append(seq_length)

        ############################################################################################
        #Somewhat tediously this block will distribute the read evenly amongst all bases in the read
        if assignment == 'all':
            if split_line[1] == '0':
                for i in range(start_loc, end_loc + 1):
                    try:
                        fwd_dicty[i] += contributions
                    except KeyError:
                        fwd_dicty[i] = contributions
            elif split_line[1] == '16':
                for i in range(start_loc, end_loc + 1):
                    try:
                        rev_dicty[i] += contributions
                    except KeyError:
                        rev_dicty[i] = contributions
        ############################################################################################
        #This block will assign the entire read to the 5' end + whatever number of integers away
        elif type(assignment) == int and assignment > 0:
            if split_line[1] == '0':
                try:
                    fwd_dicty[start_loc+assignment-1] += 1
                except KeyError:
                    fwd_dicty[start_loc+assignment-1] = 1
            elif split_line[1] == '16':
                try:
                    rev_dicty[end_loc-assignment+1] += 1
                except KeyError:
                    rev_dicty[end_loc-assignment+1] = 1
        ############################################################################################
        #Finally, this block will assign the entire read to the 3' end - whatever number of integers
        elif type(assignment) == int and assignment < 0:
            if split_line[1] == '0':
                try:
                    fwd_dicty[end_loc+assignment+1] += 1
                except KeyError:
                    fwd_dicty[end_loc+assignment+1] = 1
            elif split_line[1] == '16':
                try:
                    rev_dicty[start_loc-assignment-1] += 1
                except KeyError:
                    rev_dicty[start_loc-assignment-1] = 1
        else:
            print("Invalid location to assign reads. Please select either 'all' or positive/negative integers\
                    excluding 0 where positive denotes the number of bases from the 5' end and negative denotes\
                    the number of bases from the 3' end")
            break
    
    fwd_dicty_rpkm, rev_dicty_rpkm = normalize_dicts(fwd_dicty, rev_dicty, genome_list)

    ############################################################################################
    #Writing the basic output files
    output_file_base = sam_file_name.split('/')[-1].split('.')[0]        
    with open('../Data/{}/{}_{}_{}_{}_f.wig'.format(organism, output_file_base, length_restriction[0], length_restriction[1], assignment), 'w') as outfile:
        for pos in sorted(list(fwd_dicty.keys())):
            outfile.write('{}\t{}\n'.format(pos, fwd_dicty[pos]))
            
    with open('../Data/{}/{}_{}_{}_{}_r.wig'.format(organism, output_file_base, length_restriction[0], length_restriction[1], assignment), 'w') as outfile:
        for pos in sorted(list(rev_dicty.keys())):
            outfile.write('{}\t{}\n'.format(pos, rev_dicty[pos]))
    
    with open('../Data/{}/{}_{}_{}_{}_normed_f.wig'.format(organism, output_file_base, length_restriction[0], length_restriction[1], assignment), 'w') as outfile:
        for pos in sorted(list(fwd_dicty_rpkm.keys())):
            outfile.write('{}\t{}\n'.format(pos, fwd_dicty_rpkm[pos]))
            
    with open('../Data/{}/{}_{}_{}_{}_normed_r.wig'.format(organism, output_file_base, length_restriction[0], length_restriction[1], assignment), 'w') as outfile:
        for pos in sorted(list(rev_dicty_rpkm.keys())):
            outfile.write('{}\t{}\n'.format(pos, rev_dicty_rpkm[pos]))
    
    return mapped_len_dict

def normalize_dicts(fwd_dicty, rev_dicty, genome_list):
    total_in_CDS = []
    for genome in genome_list:
        for feature in genome.features:
            if feature.type == 'CDS':
                pos = (feature.location.start, feature.location.end)
                if feature.strand == 1:
                    for i in range(pos[0], pos[1]):
                        try:
                            total_in_CDS.append(fwd_dicty[i])
                        except KeyError:
                            total_in_CDS.append(0)
                elif feature.strand == -1:
                    for i in range(pos[0], pos[1]):
                        try:
                            total_in_CDS.append(rev_dicty[i])
                        except KeyError:
                            total_in_CDS.append(0)

    total_mapped = sum(total_in_CDS)
    fwd_dicty_rpkm = {}
    rev_dicty_rpkm = {}
    for pos in fwd_dicty:
        fwd_dicty_rpkm[pos] = (fwd_dicty[pos] * 1000000000) /total_mapped
    for pos in rev_dicty:
        rev_dicty_rpkm[pos] = (rev_dicty[pos] * 1000000000) /total_mapped
    
    return fwd_dicty_rpkm, rev_dicty_rpkm

def read_single_wiggle(wiggle_file_f, wiggle_file_r, feature_type='gene'):
    fwd_dicty = {}
    rev_dicty = {}
    with open(wiggle_file_f) as infile:
        for line in infile.readlines()[1:]:
            split_line = line.split('\t')
            fwd_dicty[int(split_line[0])] = float(split_line[1])
    with open(wiggle_file_r) as infile:
        for line in infile.readlines()[1:]:
            split_line = line.split('\t')
            rev_dicty[int(split_line[0])] = float(split_line[1])
    return fwd_dicty, rev_dicty


def read_multiple_wiggles(data_list):
    fwd_dicty_meta = {}
    rev_dicty_meta = {}
    for name, wiggle_file_f, wiggle_file_r in data_list:
        fwd_dicty, rev_dicty = read_single_wiggle(wiggle_file_f, wiggle_file_r)
        fwd_dicty_meta[name] = fwd_dicty
        rev_dicty_meta[name] = rev_dicty
    return fwd_dicty_meta, rev_dicty_meta

def get_genome_features(genome, feature_type='gene'):
    feature_dict = {}
    for feature in genome.features:
        if feature.type==feature_type:
            gene_name = feature.qualifiers['locus_tag'][0]
            feature_dict[gene_name] = feature
    return feature_dict

def get_rpkm_dict(genome_dict, fwd_dicty_meta, rev_dicty_meta):
    rpkm_dict_meta = {}
    for sample in fwd_dicty_meta.keys():
        assert sample in rev_dicty_meta.keys()
        mapped_dict = {} 
        for feature in genome_dict.values():
            pos = (feature.location.start, feature.location.end)
            sequencing = []
            if feature.strand == 1:
                for i in range(pos[0], pos[1]):
                    try:
                        sequencing.append(fwd_dicty_meta[sample][i])
                    except KeyError:
                        sequencing.append(0)
    
            elif feature.strand == -1:
                for i in range(pos[0], pos[1]):
                    try:
                        sequencing.append(rev_dicty_meta[sample][i])
                    except KeyError:
                        sequencing.append(0)
            mapped_dict[feature.qualifiers['locus_tag'][0]] = sequencing
        
        cds_sum = sum([inner for outer in list(mapped_dict.values()) for inner in outer])
        rpkm_dict_meta[sample] = {}
        for feature in mapped_dict:
            sequencing = mapped_dict[feature]
            ########################################################    
            #FPKM
            rpkm_dict_meta[sample][feature] = (sum(sequencing[:]) * 1000000000) / (len(sequencing[:]) * cds_sum)
    
    return rpkm_dict_meta

def wiggle_to_dicts(genome, sample_files, organism, file_ending, feature_type='CDS'):
    '''

    '''
    exp_dict_meta = {}
    feature_dict = {}
    sequencing_dict_meta_f = {}
    sequencing_dict_meta_r = {}
    

    for sample_file in sample_files:
        exp_dict_meta[sample_file] = {}
    
        fwd_dicty = {}
        rev_dicty = {}
        total_mapped = 0
        
        with open('../Data/{}/{}{}_f.wig'.format(organism, sample_file, file_ending)) as infile:
            for line in infile:
                split_line = line.split('\t')
                fwd_dicty[int(split_line[0])] = float(split_line[1])
        with open('../Data/{}/{}{}_r.wig'.format(organism, sample_file, file_ending)) as infile:
            for line in infile:
                split_line = line.split('\t')
                rev_dicty[int(split_line[0])] = float(split_line[1])
        
        mapped_in_features = [] 
        for feature in genome.features[:]:
            if feature.type == feature_type:
                pos = (feature.location.start, feature.location.end)
                sequencing = []
                if feature.strand == 1:
                    for i in range(pos[0], pos[1]):
                        try:
                            sequencing.append(fwd_dicty[i])
                        except KeyError:
                            sequencing.append(0)
    
                elif feature.strand == -1:
                    for i in range(pos[0], pos[1]):
                        try:
                            sequencing.append(rev_dicty[i])
                        except KeyError:
                            sequencing.append(0)
                mapped_in_features.append(sequencing)
        
        cds_sum = sum([inner for outer in mapped_in_features for inner in outer])
        
        for feature in genome.features:
            if feature.type == feature_type:
                gene_name = feature.qualifiers['locus_tag'][0]
                pos = (feature.location.start, feature.location.end)
                sequencing = []
                if feature.strand == 1:
                    for i in range(pos[0], pos[1]):
                        try:
                            sequencing.append(fwd_dicty[i])
                        except KeyError:
                            sequencing.append(0)
    
                elif feature.strand == -1:
                    for i in range(pos[0], pos[1]):
                        try:
                            sequencing.append(rev_dicty[i])
                        except KeyError:
                            sequencing.append(0)
                ########################################################    
                #FPKM
                exp_dict_meta[sample_file][gene_name] = (sum(sequencing[:]) * 1000000000) / (len(sequencing[:]) * cds_sum)
                feature_dict[gene_name] = feature

        sequencing_dict_meta_f[sample_file] = fwd_dicty
        sequencing_dict_meta_r[sample_file] = rev_dicty
        total = sum(list(fwd_dicty.values())+list(rev_dicty.values()))
        print('{}: Total reads mapped:{} Percentage mapped to feature:{}'.format(sample_file, total, (cds_sum/total)*100))
    return exp_dict_meta, feature_dict, sequencing_dict_meta_f, sequencing_dict_meta_r

def get_features_in_window(genome_features_dict, window_beg, window_end):
    feature_list = []
    for feature in genome_features_dict.values():
        if feature.location.start > window_beg and feature.location.start < window_end:
            feature_list.append(feature)
        elif feature.location.end > window_beg and feature.location.end < window_end:
            feature_list.append(feature)
        elif feature.location.start < window_beg and feature.location.end > window_end:
            feature_list.append(feature)
    return feature_list

def plot_individual_segment(window_beg, window_end, feature_list, sequencing_dict_list,\
        defined_ylimits=False, logy=False, vert_lines=[], feature_name=False, save_file=False):
    xcoords = range(window_beg, window_end)
    feature_vals_dict_plus = {}
    feature_vals_dict_minus = {}
    for label, fwd_dicty, rev_dicty, color in sequencing_dict_list:
        feature_vals_dict_plus[label] = []
        feature_vals_dict_minus[label] = []
        for position in xcoords:
            try:
                feature_vals_dict_plus[label].append(fwd_dicty[position])
            except KeyError:
                feature_vals_dict_plus[label].append(0)
            try:
                feature_vals_dict_minus[label].append(rev_dicty[position])
            except KeyError:
                feature_vals_dict_minus[label].append(0)


    fig = plt.figure(figsize=(12, 4))
    ax1 = fig.add_subplot(111)
    if feature_name:
        ax1.set_title("{}".format(feature_name), fontsize=20)
    for label, fwd_dicty, rev_dicty, color in sequencing_dict_list:
        if logy==False:
            ax1.plot(xcoords, feature_vals_dict_plus[label], label=label, c=color, linewidth=2)    
        else:
            ax1.semilogy(xcoords, feature_vals_dict_plus[label], label=label, c=color, linewidth=1)    
    for feature in feature_list:
        if feature.strand == 1:
            ax1.axvspan(feature.location.start, feature.location.end, color='k', alpha= 0.1)
    ax1.set_xlim(xcoords[0], xcoords[-1])
    ax1.legend(loc='best', fontsize=12)
    ax1.tick_params(labelsize=24)
    if logy==False:
        ax1.ticklabel_format(style='plain', useOffset=False) 
    #ax1.set_xticklabels('') 
    if defined_ylimits:
        ax1.set_ylim(0, defined_ylimits[0])
    else:
        tempy = []
        for label in feature_vals_dict_plus.keys():
            tempy.append(max(list(feature_vals_dict_plus[label])))
            tempy.append(max(list(feature_vals_dict_minus[label])))
        max_val = max(tempy)
        ax1.set_ylim(0, max_val)
    plt.xticks(rotation=45)
    plt.tight_layout()
    if logy==False:
        ax1.grid(b=True, which='both', color='0.75',linestyle='-')
    if save_file == True:
        results_path = '../Results/'+datetime.datetime.now().strftime('%Y_%m_%d/')
        if not os.path.exists(results_path):
            os.mkdir(results_path)
        plt.savefig('{}{}_avg.pdf'.format(results_path, feature_name))
        plt.close()


    #fig = plt.figure(figsize=(12, 4))
    #ax1 = fig.add_subplot(211)
    #ax2 = fig.add_subplot(212)
    #if feature_name:
    #    ax1.set_title("{}".format(feature_name), fontsize=20)
    #for label, fwd_dicty, rev_dicty in sequencing_dict_list:
    #    if logy==False:
    #        ax1.plot(xcoords, feature_vals_dict_plus[label], label=label, linewidth=1)    
    #        ax2.plot(xcoords, feature_vals_dict_minus[label], label=label, linewidth=1)    
    #    else:
    #        ax1.semilogy(xcoords, feature_vals_dict_plus[label], label=label, linewidth=1)    
    #        ax2.semilogy(xcoords, feature_vals_dict_minus[label], label=label, linewidth=1)    
    #for feature in feature_list:
    #    if feature.strand == 1:
    #        ax1.axvspan(feature.location.start, feature.location.end, color='k', alpha= 0.1)
    #    elif feature.strand == -1:
    #        ax2.axvspan(feature.location.start, feature.location.end, color='k', alpha= 0.1)
    #ax1.set_xlim(xcoords[0], xcoords[-1])
    #ax2.set_xlim(xcoords[0], xcoords[-1])
    #for vert_line in vert_lines:
    #    if vert_line[0] == 'plus':
    #       ax1.axvline(vert_line[1], c='k', linewidth=1, linestyle = '--') 
    #       ax1.axvline(vert_line[2], c='k', linewidth=1, linestyle = '--') 
    #    elif vert_line[0] == 'minus':
    #       ax2.axvline(vert_line[1], c='k', linewidth=1, linestyle = '--') 
    #       ax2.axvline(vert_line[2], c='k', linewidth=1, linestyle = '--') 
    #    else:
    #        print("mistake processing vertical lines argument")
    #ax1.legend(loc='best', fontsize=16)
    #ax1.tick_params(labelsize=12)
    #ax2.tick_params(labelsize=12)
    #if logy==False:
    #    ax1.ticklabel_format(style='plain', useOffset=False) 
    #    ax2.ticklabel_format(style='plain', useOffset=False) 
    #ax1.set_xticklabels('') 
    #if defined_ylimits:
    #    ax1.set_ylim(0, defined_ylimits[0])
    #    ax2.set_ylim(0, defined_ylimits[1])
    #else:
    #    tempy = []
    #    for label in feature_vals_dict_plus.keys():
    #        tempy.append(max(list(feature_vals_dict_plus[label])))
    #        tempy.append(max(list(feature_vals_dict_minus[label])))
    #    max_val = max(tempy)
    #    ax1.set_ylim(0, max_val)
    #    ax2.set_ylim(0, max_val)
    #plt.xticks(rotation=45)
    #ax2.invert_yaxis()
    #plt.tight_layout()
    #if logy==False:
    #    ax1.grid(b=True, which='both', color='0.75',linestyle='-')
    #    ax2.grid(b=True, which='both', color='0.75',linestyle='-')
    #if save_file == True:
    #    results_path = '../Results/'+datetime.datetime.now().strftime('%Y_%m_%d/')
    #    if not os.path.exists(results_path):
    #        os.mkdir(results_path)
    #    plt.savefig('{}{}_avg.pdf'.format(results_path, feature_name))
    #    plt.close()
    return


def length_histograms(sample_name, mapped_len_dict, organism):       
    fig = plt.figure(figsize=(10,5))
    ax1 = plt.subplot(121)
    ax1.hist(mapped_len_dict['Plus'] + mapped_len_dict['Minus'], 20)
    ax1.set_yscale('log')
    ax1.set_xlabel('length of reads', fontsize=24)
    ax1.set_ylabel('PMF', fontsize=24)
    ax1.set_title(sample_name+'\n mapped', fontsize=24)
    plt.tick_params(labelsize=18)
    ax2 = plt.subplot(122, sharey=ax1)
    ax2.hist(mapped_len_dict['Unmapped'], 20)
    ax2.set_yscale('log')
    ax2.set_xlabel('length of reads', fontsize=24)
    ax2.set_title(sample_name+'\n unmapped', fontsize=24)
    plt.tick_params(labelsize=18)
    plt.tight_layout()
    plt.savefig('../Results/'+organism+'_'+sample_name+'_lengths.pdf')
    plt.close()
    return

def analyze_rrna_percents(sample_names, generic_fq, generic_no_rrna_fq):
    '''
    The purpose of this function is to iterate through several different samples and analyze
    the percentage of reads that were mapped to rRNA. For each sample file, I need the regular .fq
    file and the rRNA depleted .fq file.
    sample_files is just a list of the sample names (i.e. ['SQ-1', 'SQ-2', etc.])
    generic_fq and generic_no_rrna_fq are just baseline file names that I 
            can replace with a given sample name (i.e. '/blahblah/SQ-X-blahblah.fq')
    '''
    import subprocess
    def file_len(fname):
        p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
                                                  stderr=subprocess.PIPE)
        result, err = p.communicate()
        if p.returncode != 0:
            raise IOError(err)
        return int(result.strip().split()[0])
    
    total_reads_dict = {}
    rrna_reads_dict = {}
    for sample_name in sample_names[:]:
        ######################################################
        ###Just counting the number of lines in the .fq files
        total_lines_fq = 0
        total_lines_no_rrna_fq = 0
        ###Open the raw .fq file, go through each line and enumerate
        file_name = generic_fq.replace('XXXX', sample_name)
        total_reads = file_len(file_name)/4
        file_name = generic_no_rrna_fq.replace('XXXX', sample_name)
        sans_rrna_reads = file_len(file_name)/4
        rrna_reads = total_reads - sans_rrna_reads
        total_reads_dict[sample_name] = total_reads
        rrna_reads_dict[sample_name] = rrna_reads
        
    return total_reads_dict, rrna_reads_dict



def wiggle_to_features(fwd_wiggle, rev_wiggle, genome, sample_name, length_restriction=200, coverage_restriction=10, utr_length=50):
    feature_sequencing = {}
    fwd_dicty = {}
    rev_dicty = {}
    total_mapped = 0 
    fwd_dicty, rev_dicty = read_single_wiggle(fwd_wiggle, rev_wiggle)
    
    for feature in genome.features[:]:
        if feature.type == 'CDS':
            gene_name = feature.qualifiers['locus_tag'][0]
            pos = (feature.location.start-utr_length, feature.location.end+utr_length) ####NOTE: including UTRs here!!!!
            sequencing = []
            if feature.strand == 1:
                for i in range(pos[0], pos[1]):
                    try:
                        sequencing.append(fwd_dicty[i])
                    except KeyError:
                        sequencing.append(0)
                feature_sequencing[gene_name] = sequencing
            elif feature.strand == -1: 
                for i in range(pos[0], pos[1]):
                    try:
                        sequencing.append(rev_dicty[i])
                    except KeyError:
                        sequencing.append(0)
                feature_sequencing[gene_name] = sequencing[::-1]
    
    print('{} initially had {} genes'.format(sample_name, len(feature_sequencing)))
    to_delete = []
    for gene_name in feature_sequencing:
        if len(feature_sequencing[gene_name]) < length_restriction:
            to_delete.append(gene_name)
        elif np.percentile(feature_sequencing[gene_name], 100-coverage_restriction) <= 0: ###Important parameter here, genes to discard based on coverage
            to_delete.append(gene_name)
    for gene_name in to_delete:
        del feature_sequencing[gene_name]
    print('... it now has {} after coverage and length pruning'.format(len(feature_sequencing)))
    return feature_sequencing 
    
def meta_gene_analysis(feature_sequencing, sample_name, winsorize=True, save_file_name=False):
    meta = []
    x_vals = []
    for i in range(150):
        tempy = []
        for gene, reads in feature_sequencing.items():
            tempy.append(reads[i])
            
        ######IMPORTANT PARAMETER HERE
        if winsorize:
            meta.append(np.mean(stats.mstats.winsorize(tempy, limits=0.1))) ###Winsorize
        else:
            meta.append(np.mean(tempy)) ###Don't Winsorize
        x_vals.append(i-50)

    fig = plt.figure(figsize=(12,4))
    ax1 = fig.add_subplot(121)
    ax1.plot(x_vals, meta)
    ax1.axvline(0, c='k', linestyle='--')
    ax1.set_title("{}: 5' alignment\n n={} genes".format(sample_name, len(feature_sequencing)), fontsize=24)
    ax1.tick_params(labelsize=16)
    ax1.grid(b=True, which='both', color='0.75',linestyle='-')
    ax1.set_xlabel('Location relative to \nstart codon (nt)', fontsize=18)
    meta = []
    x_vals = []
    for i in range(150, 0, -1):
        tempy = []
        for gene, reads in feature_sequencing.items():
            tempy.append(reads[-1*i])
        if winsorize:
            meta.append(np.mean(stats.mstats.winsorize(tempy, limits=0.1))) ###Winsorize
        else:
            meta.append(np.mean(tempy)) ###Don't Winsorize
        x_vals.append((-1*i)+50)

    ax2 = fig.add_subplot(122, sharey=ax1)
    ax2.plot(x_vals, meta)
    ax2.axvline(0, c='k', linestyle='--')
    ax2.set_title("{}: 3' alignment\n n={} genes".format(sample_name, len(feature_sequencing)), fontsize=24)
    ax2.tick_params(labelsize=16)
    #ax2.set_yticklabels('')
    ax2.grid(b=True, which='both', color='0.75',linestyle='-')
    ax2.set_xlabel('Location relative to \nstop codon (nt)', fontsize=18)
    fig.tight_layout()
    
    if save_file_name:
        results_path = '../Results/'+datetime.datetime.now().strftime('%Y_%m_%d/')
        if not os.path.exists(results_path):
            os.mkdir(results_path)
        plt.savefig('{}{}'.format(results_path, save_file_name))
        
    return





def plot_correlations_single(first_data_label, second_data_label, df, save_file_name=False):
    plt.figure()
    plt.title("Spearman's rho = {}".format(stats.spearmanr(df[first_data_label], df[second_data_label])[0]))
    plt.loglog(df[first_data_label], df[second_data_label], 'bo')
    plt.xlabel('{}'.format(first_data_label))
    plt.ylabel('{}'.format(second_data_label))
    if save_file_name:
        results_path = '../Results/'+datetime.datetime.now().strftime('%Y_%m_%d/')
        if not os.path.exists(results_path):
            os.mkdir(results_path)
        plt.savefig('{}{}'.format(results_path, save_file_name))
    return

def plot_correlations_combined(first_data_label, second_data_label, df, p_val_threshold = 0.001, save_file_name=False):
    colors = []
    for i in df['T-test(p-value)']:
        if i < p_val_threshold:
            colors.append('r')
        else:
            colors.append('k')
    
    max_vals = max([max(df[first_data_label]), max(df[second_data_label])])  
    fig = plt.figure(figsize=(12,8))
    ax1 = fig.add_subplot(111)
    ax1.scatter(df[first_data_label], df[second_data_label], c=colors, alpha=0.5);
    ax1.set_title("Spearman's rho = {}".format(stats.spearmanr(df[first_data_label], df[second_data_label])[0]))
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_xlim(10e-3, max_vals*10)
    ax1.set_ylim(10e-3, max_vals*10)
    ax1.tick_params(labelsize=20)
    ax1.set_xlabel('RPKM {}'.format(first_data_label), fontsize=24)
    ax1.set_ylabel('RPKM {}'.format(second_data_label), fontsize=24)
    if save_file_name:
        results_path = '../Results/'+datetime.datetime.now().strftime('%Y_%m_%d/')
        if not os.path.exists(results_path):
            os.mkdir(results_path)
        plt.savefig('{}{}'.format(results_path, save_file_name))

    fig = plt.figure(figsize=(12,8))
    ax1 = fig.add_subplot(111)
    ax1.scatter(df['log2(fold change)'], -1*np.log10(df['T-test(p-value)']), c=colors, alpha=0.3)
    ax1.axvline(1, color='c', linestyle='--')
    ax1.axvline(-1, color='c', linestyle='--')
    ax1.set_xlabel('log2(fold change)', fontsize=24)
    ax1.set_ylabel('-log10(p-value)', fontsize=24)
    ax1.tick_params(labelsize=20)
    if save_file_name:
        results_path = '../Results/'+datetime.datetime.now().strftime('%Y_%m_%d/')
        if not os.path.exists(results_path):
            os.mkdir(results_path)
        plt.savefig('{}Volcano_{}'.format(results_path, save_file_name))
    return

def plot_pvals(window_beg, window_end, feature_list, control_labels, treatment_labels, sequencing_dict_meta_f_normed, sequencing_dict_meta_r_normed, window_size=100, step_size=10):
    ctrl_fpkm_avg_plus = []
    ctrl_fpkm_avg_minus = []
    treatment_fpkm_avg_plus = []
    treatment_fpkm_avg_minus = []
    
    p_vals_plus = []
    p_vals_minus = []
    
    xcoords = []
    for i in range(window_beg, window_end+step_size, step_size):
        treatment_vals_plus = []
        control_vals_plus = []
        treatment_vals_minus = []
        control_vals_minus = []
        for label in sequencing_dict_meta_f_normed.keys():
            temp_plus = []
            temp_minus = []
            for j in range(i-int((0.5*window_size)), i+int(0.5*window_size)):
                try:
                    temp_plus.append(sequencing_dict_meta_f_normed[label][j])
                except KeyError:
                    temp_plus.append(0)
                try:
                    temp_minus.append(sequencing_dict_meta_r_normed[label][j])
                except KeyError:
                    temp_minus.append(0)
                            
            fpkm_plus = np.mean(temp_plus)
            fpkm_minus = np.mean(temp_minus)
    
            if label in control_labels:
                control_vals_plus.append(fpkm_plus)
                control_vals_minus.append(fpkm_minus)
    
            elif label in treatment_labels:
                treatment_vals_plus.append(fpkm_plus)
                treatment_vals_minus.append(fpkm_minus)
        
        ctrl_fpkm_avg_plus.append(np.mean(control_vals_plus))
        ctrl_fpkm_avg_minus.append(np.mean(control_vals_minus))
        treatment_fpkm_avg_plus.append(np.mean(treatment_vals_plus))
        treatment_fpkm_avg_minus.append(np.mean(treatment_vals_minus))
        t, p = stats.ttest_ind(control_vals_plus, treatment_vals_plus, equal_var=False)
        p_vals_plus.append(p)
        t, p = stats.ttest_ind(control_vals_minus, treatment_vals_minus, equal_var=False)
        p_vals_minus.append(p)
        xcoords.append(i)
    
    fig = plt.figure(figsize=(16.1,4))
    ax1 = fig.add_subplot(211)
    ax1.semilogy(xcoords, p_vals_plus, color='k')
    ax1.axhline(0.01, color='r')
    ax1.axhline(0.001, color='r')
    for feature in feature_list:
        if feature.strand == 1:
            ax1.axvspan(feature.location.start, feature.location.end, color='k', alpha= 0.1)
    ax1.set_xlim(xcoords[0], xcoords[-1])
    #ax1.ticklabel_format(style='plain', useOffset=False) 
    ax1.tick_params(labelsize=16)
    ax1.set_xticklabels('')
    ax2 = fig.add_subplot(212)
    ax2.semilogy(xcoords, p_vals_minus, color='k')
    for feature in feature_list:
        if feature.strand == -1:
            ax2.axvspan(feature.location.start, feature.location.end, color='k', alpha= 0.1)
    ax2.axhline(0.01, color='r')
    ax2.axhline(0.001, color='r')
    ax2.set_xlim(xcoords[0], xcoords[-1])
    #ax2.ticklabel_format(style='plain', useOffset=False) 
    ax2.tick_params(labelsize=16)
    plt.xticks(rotation=45)

    return
