#!/usr/bin/python

import numpy as np, numpy.linalg
import argparse, math, sys
from src import io

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
    snp_idx,all_snp_leg = io.load_legend(leg_file)

    # load zscore file
    snp_beta,all_snp_beta = io.load_beta(zsc_file)
    
    # load partition
    part = io.load_partition(part_file)
    
    # output eigenvalue and projection squared
    #output_eig_prjsq(chrom, snp_idx, all_snp_leg, snp_beta,
        #all_snp_beta, part, ref_file, out_file)

# output eigenvalue decomposition and projection squared
def output_eig_prjsq(chrom, snp_idx, all_snp_leg, snp_beta,
    all_snp_beta, part, ref_file, out_file):
    
    # open files to read and write
    out_file_info = out_file+'_chr'+chrom+'.info'
    out_file_proj = out_file+'_chr'+chrom+'.prjsq'
    out_file_eig = out_file+'_chr'+chrom+'.eig'
    out_file_info = open(out_file_info, 'w')
    out_file_proj = open(out_file_proj, 'w')
    out_file_eig = open(out_file_eig, 'w')
    ref_file = open(ref_file)
    
    # iterate through locus
    line_idx = 0
    start_idx = 0
    for locus in part:

        # find snps in a locus
        locus_snp,end_idx = get_locus_snp(all_snp_beta, locus, start_idx)
        start_idx = end_idx

        # get beta for snps in the locus
        locus_beta = get_locus_snp_beta(snp_beta, locus_snp)
        
        # load reference panel for the locus
        load_line_idx = get_load_line_idx(snp_idx, locus_snp)
        snp_hap,line_idx,last_fp = load_ref_file(ref_file, load_line_idx,
            all_snp_leg, line_idx)
        gens = get_locus_refpanel(snp_hap, locus_snp)
        
        # check for empty locus
        if(len(locus_beta) == 0):
            out_file_info.write('%d\t%d\t%d\t%d\t%.1f\n'
                % (locus[0], locus[1], 0, 0, 0.0))
            out_file_eig.write('%.8f\n' % (0.0))
            out_file_proj.write('%.8f\n' % (0.0))
            continue
        
        # compute window ld and its decomposition
        locus_ld = get_ld(gens)
        locus_ld_rank = np.linalg.matrix_rank(locus_ld)
        ld_w,ld_v = eig_decomp(locus_ld)

        # write window info
        locus_n = [elem[2] for elem in locus_snp]
        locus_info = '%d\t%d\t%d\t%d\t%.1f\n' % (locus[0], locus[1],
            gens.shape[0], locus_ld_rank, np.mean(np.array(locus_n)))
        out_file_info.write(locus_info)

        # write eigen values and squared projection
        locus_eig = ''
        locus_proj = ''
        for i in xrange(locus_ld_rank):
            locus_eig += '%.8f\t' % (ld_w[i,0])
            beta_vec = np.matrix(locus_beta).T
            eig_vec = np.matrix(np.real(ld_v[:,i]))
            locus_proj += '%.8f\t' % (((beta_vec.T*eig_vec)[0,0])**2.0)
        out_file_eig.write(locus_eig+'\n')
        out_file_proj.write(locus_proj+'\n')

    # close files
    ref_file.close()
    out_file_info.close()
    out_file_eig.close()
    out_file_proj.close()

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
