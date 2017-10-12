#!/usr/bin/python
# (c) 2017-2022 Huwenbo Shi

import numpy as np
import pandas as pd
import argparse, math, sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec

# main function
def main():
    
    # get command line input
    args = get_command_line()

    # check if number of files matches number of trait names
    if len(args.local_hsq_files) != len(args.trait_names):
        sys.stderr.write('Error: Number of local SNP-heritability estimate'\
            'files does not match number of trait names\n')
        sys.exit()
    ntraits = len(args.trait_names)

    # load results
    all_local_hsq = []
    for i in xrange(ntraits):
        tmp = pd.read_table(args.local_hsq_files[i], sep='\t')
        if args.no_negative:
            tmp.loc[tmp['local_h2g']<0, 'local_h2g'] = 0.0
        all_local_hsq.append(tmp)
        
    # sort local SNP-heritability by descending order and get total hsq and
    # total number of snps
    all_total_hsq = []
    all_total_nsnp = []
    all_local_hsq_sorted = []
    all_nsnp_sorted = []
    for i in xrange(ntraits):
        
        # sort local SNP-heritability by descending order
        idx = np.argsort(-all_local_hsq[i]['local_h2g'])
        all_local_hsq_sorted.append(all_local_hsq[i]['local_h2g'][idx]\
            .reset_index(drop=True))
        all_nsnp_sorted.append(all_local_hsq[i]['num_snp'][idx]\
            .reset_index(drop=True))
        
        # get total SNP-heritability and total number of SNPs
        all_total_hsq.append(np.sum(all_local_hsq[i]['local_h2g']))
        all_total_nsnp.append(np.float(np.sum(all_local_hsq[i]['num_snp'])))
    
    # get the proportions
    all_xval = [np.array([]) for i in xrange(ntraits)]
    all_yval = [np.array([]) for i in xrange(ntraits)]
    all_yerr = [np.array([]) for i in xrange(ntraits)]

    # iterate through the traits
    nloci = all_local_hsq_sorted[0].shape[0]
    for i in xrange(ntraits):
        
        # iterate through the loci
        for j in xrange(nloci):
            hsq_sum = np.sum(all_local_hsq_sorted[i][0:j+1])
            hsq_frac = hsq_sum / all_total_hsq[i]
            snp_frac = np.sum(all_nsnp_sorted[i][0:j+1]) / all_total_nsnp[i]
            all_xval[i] = np.append(all_xval[i], snp_frac)
            all_yval[i] = np.append(all_yval[i], hsq_frac)

            # get jack knife standard error
            if args.show_se:
                jk_est = []
                for k in xrange(nloci):
                    total_hsq_jk = all_total_hsq[i]-all_local_hsq_sorted[i][k]
                    hsq_sum_jk = hsq_sum
                    if k <= j:
                        hsq_sum_jk -= all_local_hsq_sorted[i][k]
                    jk_est.append(hsq_sum_jk/total_hsq_jk)
                jk_est = np.array(jk_est)
                se = np.sqrt((nloci-1)*np.mean(np.square(jk_est-hsq_frac)))
                all_yerr[i] = np.append(all_yerr[i], se)

    # plot the ratio vs ratio plots
    fig, ax = plt.subplots()
    color_cycle = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
        "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
    for i in xrange(ntraits):
        ax.plot(all_xval[i], all_yval[i], color=color_cycle[i%ntraits],
            label=args.trait_names[i])
        
        if args.show_se:
            ax.fill_between(all_xval[i], all_yval[i]-1.96*all_yerr[i],
                all_yval[i]+1.96*all_yerr[i], alpha=0.2,
                edgecolor=color_cycle[i%ntraits],
                facecolor=color_cycle[i%ntraits], linewidth=0.0)

    # plot the y=x line
    ax.plot([0, 1.01], [0, 1.01], 'k--', label='y = x')

    # set legend
    plt.legend(loc=4)

    # set labels
    ax.set_xlabel('Percent of genome covered')
    ax.set_ylabel('Percent of total SNP-heritability')
    
    # clip figure
    ax.set_xlim([0, 1.00])
    ax.set_ylim([0, 1.01])

    # save figure
    plt.savefig(args.out_file, bbox_inches='tight')

# get command line
def get_command_line():
    
    parser = argparse.ArgumentParser(description='Contrast polygenicity' \
        'between multiple traits using estimates of local SNP-heritability')

    parser.add_argument('--local-hsqg-est', dest='local_hsq_files', type=str,
        help='A list of local files containing local SNP-heritability '\
        'estimates', nargs='+', required=True)
    
    parser.add_argument('--trait-names', dest='trait_names', type=str,
        help='A list of trait names', nargs='+', required=True)

    parser.add_argument('--no-negative', dest='no_negative',
        help='Set the negative local SNP-heritability estimates to be zero',
        default=False, action='store_true', required=False)

    parser.add_argument('--show-se', dest='show_se',
        help='Show standard error band on the plots', default=False,
        action='store_true', required=False)

    parser.add_argument('--out', dest='out_file', type=str,
        help='Output figure file name', required=True)
    
    args = parser.parse_args()
    
    return args

if(__name__ == '__main__'):
    main()
