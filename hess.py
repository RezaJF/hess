#!/usr/bin/python
# (c) 2016-2021 Huwenbo Shi

import numpy as np, numpy.linalg
import argparse, math, sys, logging, time
from src import io, step1, step2

# main function
def main():

    # get command line
    args = get_command_line()

    # get step 1 arguments
    zsc_file = args.zscore_file
    leg_file = args.legend_file
    out_file_step1 = args.output_file
    ref_file = args.reference_panel
    part_file = args.partition_file
    chrom = args.chrom

    # get step 2 arguments
    prefix = args.prefix
    num_eig = args.k
    out_file_step2 = args.out
    gc = args.lambda_gc
    if(args.tot_h2g is not None):
        tot_h2g = args.tot_h2g[0]
        tot_h2g_se = args.tot_h2g[1]
    sense_thres_joint = args.sense_threshold_joint
    sense_thres_indep = args.sense_threshold_indep
    eig_thres = args.eig_threshold   
 
    ##########     run step 1     ##########
    if(zsc_file        is not None and 
       leg_file        is not None and
       out_file_step1  is not None and
       ref_file        is not None and
       part_file       is not None and
       chrom           is not None and
       prefix          is     None and
       out_file_step2  is     None and
       args.tot_h2g    is     None):

        # starting the log
        logging.basicConfig(filename=out_file_step1+'_chr'+chrom+'.log',
            level=logging.INFO, format='%(message)s', filemode="w")
        cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        logging.info('Command started at: %s' % cur_time)
        logging.info('Command issued: %s' % ' '.join(sys.argv))

        # load snp in legend
        refpanel_snp_idx, refpanel_leg = io.load_legend(leg_file)
        logging.info('Number of SNPs in reference panel: '
                     '%d' % len(refpanel_leg))
        
        # load zscore file
        snp_beta,snp_beta_info = io.load_beta(zsc_file)
        logging.info('Number of SNPs in Z-score file: '
                     '%d' % len(snp_beta))
        
        # load partition
        part = io.load_partition(part_file)
        logging.info('Number of loci in partition file: '
                     '%d' % len(part))
        
        # output eigenvalue and projection squared
        step1.output_eig_prjsq(chrom, refpanel_snp_idx, refpanel_leg,
            snp_beta, snp_beta_info, part, ref_file, out_file_step1)
        
        # ending the log
        cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
        logging.info('Command finished at: %s' % cur_time)

    ##########     run step 2     ##########
    elif(zsc_file        is     None and 
         leg_file        is     None and
         out_file_step1  is     None and
         ref_file        is     None and
         part_file       is     None and
         chrom           is     None and
         prefix          is not None and
         out_file_step2  is not None):
        
        # load step1
        locus_info,all_eig,all_prj = io.load_step1(prefix)
        
        # estimate h2g jointly when total h2g is not provided
        if(args.tot_h2g is None):
            all_h2g,raw_est = step2.get_local_h2g_joint(locus_info, all_eig,
                                        all_prj, num_eig, eig_thres,
                                        sense_thres_joint, gc)
            all_var = step2.get_var_est_joint(locus_info, all_h2g)
        # estimate h2g independently when total h2g is provided
        else:
            all_h2g,raw_est = step2.get_local_h2g_indep(locus_info, all_eig,
                                        all_prj, num_eig, eig_thres,
                                        sense_thres_indep, gc, tot_h2g)
            all_var = step2.get_var_est_indep(locus_info, all_h2g, tot_h2g_se)

        # write output
        step2.output_local_h2g(out_file_step2, locus_info, raw_est,
                all_h2g, all_var)
    
    ##########     unrecognized     ##########
    else:
        print 'unrecognized option'
        sys.exit(1)

# get command line
def get_command_line():
    parser = argparse.ArgumentParser(description='Compute eigenvalues \
                    of LD matrices, and squared projections of effect \
                    size vector onto corresponding eigenvectors')

    # step 1 arguments
    parser.add_argument('--zscore-file', dest='zscore_file', type=str,
                   help='Z-score file', required=False)
    parser.add_argument('--chrom', dest='chrom', type=str,
                   help='Chromosome number', required=False)
    parser.add_argument('--output-file', dest='output_file', type=str,
                   help='Output file name', required=False)
    parser.add_argument('--reference-panel', dest='reference_panel', type=str,
                   default=None, help='Reference panel file', required=False)
    parser.add_argument('--legend-file', dest='legend_file', type=str,
                   help='Legend file', required=False)
    parser.add_argument('--partition-file', dest='partition_file', type=str,
                   help='Partition file', required=False)
    
    
    # step 2 arguments
    parser.add_argument('--prefix', dest='prefix', type=str,
                   help='Prefix used for step 1', required=False)
    parser.add_argument('--out', dest='out', type=str,
                   help='Output file name', required=False)
    parser.add_argument('--k', dest='k', type=int, default=50,
                   help='Maximum number of eigenvectors to use (default 50)')
    parser.add_argument('--lambda_gc', dest='lambda_gc', type=float,
                   default=None, help='Genomic control factor (will be \
                   estimated if not provided)')
    parser.add_argument('--tot-h2g', dest='tot_h2g', nargs=2, type=float,
                   help='Total trait SNP heritability', required=False)
    parser.add_argument('--sense-threshold-joint',
                   dest='sense_threshold_joint', type=float,
                   default=2.0, help='Sensitivity threshold on \
                   total h2g estimates, used when tot_h2g is not provided \
                   (default 2.0)')
    parser.add_argument('--sense-threshold-indep',
                   dest='sense_threshold_indep', type=float,
                   default=0.5, help='Sensitivity threshold on \
                   total h2g estimates, used when tot_h2g is provided \
                   (default 0.5)')
    parser.add_argument('--eig-threshold', dest='eig_threshold', type=float,
                   default=1.0, help='Eigenvalue threshold (default 1.0)') 
 
    args = parser.parse_args()
    
    return args

if(__name__ == '__main__'):
    main()
