import sys, math
import numpy as np

"""
description:
    load legend file for reference panel
return:
    1. a dictionary that maps snp id with line number
    2. a list of snp ids in the order they appear in the legend
arguments:
    1. legend_file_name (str) - path to the legend file
""" 
def load_legend(legend_file_name):
    all_snps = []
    legend_file = open(legend_file_name, 'r')
    snp_idx = dict()
    idx = 0
    for line in legend_file:
        cols = line.strip().split()
        snpid = cols[0]
        all_snps.append(snpid)
        snp_idx[snpid] = idx
        idx += 1
    legend_file.close()
    return (snp_idx, all_snps)


"""
description:
    load effect size (beta) from z-score file
return:
    1. a dictionary that maps snp id with beta
    2. a list of (snp id, position, sample size) in the order they
       appear in the z-score file
arguments:
    1. zscore_file_name (str) - path to the z-score file
"""
def load_beta(zscore_file_name):
    all_snps = []
    snp_beta = dict()
    zscore_file = open(zscore_file_name, 'r')
    for line in zscore_file:
        cols = line.strip().split()
        pos = int(cols[1])
        zscore = float(cols[4])
        n = float(cols[5])
        beta = zscore / math.sqrt(n)
        snp_beta[cols[0]] = beta
        all_snps.append((cols[0], pos, n))
    zscore_file.close()
    return (snp_beta, all_snps)


"""
description:
    load partition file
return:
    1. a list of (start position, end position)
arguments:
    1. partition_file_name (str) - path to the partition file 
"""
def load_partition(partition_file_name):
    partition = []
    first_line_read = False
    partition_file = open(partition_file_name, 'r')
    for line in partition_file:
        # skip first line
        if(first_line_read == False):
            first_line_read = True
            continue
        cols = line.strip().split()
        start_pos = int(cols[1])
        end_pos = int(cols[2])-1
        partition.append((start_pos, end_pos))
    partition_file.close()
    return partition


"""
description:
    load specific lines in the reference panel
return:
arguments:
"""
def load_reference_panel(ref_panel_file, lines_to_load,
                         legend, start_line):
    snp_idx = dict()
    ref_data = np.array([])
    end_line = start_line
    num_snp_to_load = len(lines_to_load)
    num_snp_loaded = 0
    idx = 0
    while(num_snp_loaded < num_snp_to_load):
        line = ref_panel_file.readline()
        if(not line): break
        if(end_line in lines_to_load):
            snp = legend[new_start_line]
            cols = line.strip().split()
            ref_data = np.append(ref_data, np.array(cols).astype(float))
            snp_idx[snp] = idx
            num_snp_loaded += 1
            idx += 1
        end_line += 1
    return (ref_data, snp_idx, end_line)
