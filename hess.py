#!/usr/bin/python

import numpy as np, numpy.linalg
import argparse, math, sys
from src import io
from src import step1

# main function
def main():

    # get command line
    args = get_command_line()  
    zsc_file = args.zscore_file
    leg_file = args.legend_file
    out_file = args.output_file
    ref_file = args.reference_panel
    part_file = args.partition_file
    chrom = args.chrom
    
    # load snp in legend
    refpanel_snp_idx, refpanel_leg = io.load_legend(leg_file)

    # load zscore file
    snp_beta,snp_beta_info = io.load_beta(zsc_file)
    
    # load partition
    part = io.load_partition(part_file)
    
    # output eigenvalue and projection squared
    step1.output_eig_prjsq(chrom, refpanel_snp_idx, refpanel_leg, snp_beta,
        snp_beta_info, part, ref_file, out_file)

# get command line
def get_command_line():
    parser = argparse.ArgumentParser(description='Compute eigenvalues \
                    of LD matrices, and squared projections of effect \
                    size vector onto corresponding eigenvectors')
    parser.add_argument('--zscore-file', dest='zscore_file', type=str,
                   help='Z-score file', required=True)
    parser.add_argument('--chrom', dest='chrom', type=str,
                   help='Chromosome number', required=True)
    parser.add_argument('--output-file', dest='output_file', type=str,
                   help='Output file name', required=True)
    parser.add_argument('--reference-panel', dest='reference_panel', type=str,
                   default=None, help='Reference panel file', required=True)
    parser.add_argument('--legend-file', dest='legend_file', type=str,
                   help='Legend file', required=True)
    parser.add_argument('--partition-file', dest='partition_file', type=str,
                   help='Partition file', required=True)
    args = parser.parse_args()
    return args

if(__name__ == '__main__'):
    main()
