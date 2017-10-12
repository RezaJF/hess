#!/usr/bin/python
# (c) 2016-2021 Huwenbo Shi

import numpy as np
import pandas as pd
import argparse, math, sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib import gridspec
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['ytick.labelsize'] = 14
matplotlib.rc('font', family='sans-serif') 

alpha = 0.05

def main():

    # get command line input
    args = get_command_line()
    hsq = pd.read_table(args.local_hsqg_file)
    trait = args.trait_name.upper()
    num_loci = hsq.shape[0]
    thres = alpha / num_loci

    # create indices
    index = hsq.index.values
    even_chr_idx = index[np.where(hsq['chr']%2 == 0)]
    odd_chr_idx = index[np.where(hsq['chr']%2 == 1)]
    sig_even_chr_idx = index[np.where((hsq['chr']%2 == 0)&(hsq['p']<thres))]
    sig_odd_chr_idx = index[np.where((hsq['chr']%2 == 1)&(hsq['p']<thres))]
    
    # set up canvas
    fig, ax = plt.subplots(figsize=(20, 3))
    fig.tight_layout()

    # make bar plot for hsq #################################################
    ax.bar(even_chr_idx, hsq['local_h2g'][even_chr_idx], 1,
        color='whitesmoke', edgecolor='whitesmoke', linewidth=0.1)
    ax.bar(odd_chr_idx, hsq['local_h2g'][odd_chr_idx], 1,
        color='gainsboro', edgecolor='gainsboro', linewidth=0.1)
    
    # plot significant loci
    ax.bar(sig_even_chr_idx, hsq['local_h2g'][sig_even_chr_idx], 1,
        color='r', edgecolor='r', linewidth=0.1)
    ax.bar(sig_odd_chr_idx, hsq['local_h2g'][sig_odd_chr_idx], 1,
        color='b', edgecolor='b', linewidth=0.1)
    
    # get tick location
    tick_index = [0]*22
    for i in xrange(22):
        tick_index[i] = int(np.mean(index[np.where(hsq['chr']==(i+1))[0]]))

    # modify ticks
    ax.set_xticks(tick_index)
    ax.set_xticklabels(range(1,num_loci+1))
    ax.tick_params(axis=u'both', which=u'both',length=0)

    # tuning figure
    ax.set_xlim(0, num_loci)
    ax.set_ylim(0, 0.002)
    ax.set_ylabel(trait, fontsize=20)
    ax.set_title("local SNP-heritability", fontsize=20)
    tick_int = 0.002/5.0
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_int))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%g'))
    # make bar plot for hsq #################################################


    # save figure ###########################################################
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%g'))
    plt.savefig(args.out, bbox_inches='tight')

# get command line
def get_command_line():
    
    parser = argparse.ArgumentParser(description='Make manhattan plot for '\
        'local SNP-heritability estimates')
    
    parser.add_argument('--local-hsqg-est', dest='local_hsqg_file', type=str,
        help='Local SNP-heritability estimation results', required=True)
    
    parser.add_argument('--trait-name', dest='trait_name', type=str,
        help='Name of the trait', required=True)
    
    parser.add_argument('--out', type=str, dest='out', required=True,
        help='Output file name')

    args = parser.parse_args()
    
    return args

if(__name__ == '__main__'):
    main()
