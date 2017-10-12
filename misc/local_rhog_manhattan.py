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
    gcov = pd.read_table(args.local_rhog_file)
    hsq1 = pd.read_table(args.local_hsqg_file[0])
    hsq2 = pd.read_table(args.local_hsqg_file[1])
    trait1, trait2 = args.trait_names
    trait1 = trait1.upper()
    trait2 = trait2.upper()
    num_loci = hsq1.shape[0]
    thres = alpha / num_loci

    # get the correlation estimates
    gcor = np.array([0.0]*num_loci)
    pos_idx = (hsq1['local_h2g']>0) & (hsq2['local_h2g']>0)
    pos_idx = np.where(pos_idx == True)[0]
    denom = np.multiply(hsq1['local_h2g'][pos_idx], hsq2['local_h2g'][pos_idx])
    denom = np.sqrt(denom)
    gcor[pos_idx] = np.divide(gcov['local_rhog'][pos_idx], denom)
   
    # create indices
    index = gcov.index.values
    even_chr_idx = index[np.where(gcov['chr']%2 == 0)]
    odd_chr_idx = index[np.where(gcov['chr']%2 == 1)]
    sig_even_chr_idx = index[np.where((gcov['chr']%2 == 0)&(gcov['p']<thres))]
    sig_odd_chr_idx = index[np.where((gcov['chr']%2 == 1)&(gcov['p']<thres))]
    
    # set up canvas
    fig = plt.figure(figsize=(20, 6))

    gs = gridspec.GridSpec(4, 1, height_ratios=[2, 4, 1.0, 1.0]) 
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])
    ax3 = plt.subplot(gs[3])
    
    fig.tight_layout()

    # make bar plot for gcov #################################################
    ax1.bar(even_chr_idx, gcov['local_rhog'][even_chr_idx], 1,
        color='whitesmoke', edgecolor='whitesmoke', linewidth=0.1)
    ax1.bar(odd_chr_idx, gcov['local_rhog'][odd_chr_idx], 1,
        color='gainsboro', edgecolor='gainsboro', linewidth=0.1)
    
    # plot significant loci
    ax1.bar(sig_even_chr_idx, gcov['local_rhog'][sig_even_chr_idx], 1,
        color='r', edgecolor='r', linewidth=0.1)
    ax1.bar(sig_odd_chr_idx, gcov['local_rhog'][sig_odd_chr_idx], 1,
        color='b', edgecolor='b', linewidth=0.1)
    
    # get tick location
    tick_index = [0]*22
    for i in xrange(22):
        tick_index[i] = int(np.mean(index[np.where(gcov['chr']==(i+1))[0]]))

    # modify ticks
    ax1.set_xticks(tick_index)
    ax1.set_xticklabels(range(1,num_loci+1))
    ax1.tick_params(axis=u'both', which=u'both',length=0)

    # tuning figure
    ax1.set_xlim(0, num_loci)
    ax1.set_ylabel(trait1+' & '+trait2, fontsize=18)
    ax1.yaxis.set_label_coords(-0.04,1.02)
    ax1.set_title("local genetic covariance", fontsize=18)
    tick_int = (np.max(gcov['local_rhog'])-np.min(gcov['local_rhog']))/5.0
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(tick_int))
    ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.4f'))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0)) 
    # make bar plot for gcov #################################################


    # make bar plot for gcor ################################################
    ax0.bar(even_chr_idx, gcor[even_chr_idx], 1, color='whitesmoke',
        edgecolor='whitesmoke', linewidth=0.1)
    ax0.bar(odd_chr_idx, gcor[odd_chr_idx], 1, color='gainsboro',
        edgecolor='gainsboro', linewidth=0.1)
    
    # plot significant loci
    ax0.bar(sig_even_chr_idx, gcor[sig_even_chr_idx], 1, color='r',
        edgecolor='r', linewidth=0.1)
    ax0.bar(sig_odd_chr_idx, gcor[sig_odd_chr_idx], 1, color='b',
        edgecolor='b', linewidth=0.1)

    # modify ticks
    ax0.set_xticks(tick_index)
    ax0.set_xticklabels(range(1,num_loci+1))
    ax0.tick_params(axis=u'both', which=u'both',length=0)

    # tuning figure
    ax0.set_xlim(0, num_loci)
    ax0.set_ylim(-1, 1)
    ax0.set_title("local genetic correlation", fontsize=18)
    tick_int = 3.0/3.0
    ax0.yaxis.set_major_locator(ticker.MultipleLocator(tick_int))
    ax0.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    # make bar plot for gcor #################################################


    # make bar plot for h2g ##################################################
    ax2.bar(odd_chr_idx, hsq1['local_h2g'][odd_chr_idx], 1, color='b',
        edgecolor='b', linewidth=0.1)
    ax2.bar(even_chr_idx, hsq1['local_h2g'][even_chr_idx], 1, color='r',
        edgecolor='r', linewidth=0.1)

    # modify ticks
    ax2.set_xticks(tick_index)
    ax2.set_xticklabels(range(1,num_loci+1))
    ax2.tick_params(axis=u'both', which=u'both',length=0)

    # tuning figure
    ax2.set_xlim(0, num_loci)
    ax2.set_ylabel(trait1, fontsize=18)
    ax2.set_ylim(0.0, 0.002)
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.001))
    ax2.set_title("local SNP-heritability", fontsize=18)
    plt.ticklabel_format(style='plain', axis='y', scilimits=(0,0))
    # make bar plot for h2g ##################################################


    # make bar plot for h2g ##################################################
    ax3.bar(odd_chr_idx, hsq2['local_h2g'][odd_chr_idx], 1, color='b',
        edgecolor='b', linewidth=0.1)
    ax3.bar(even_chr_idx, hsq2['local_h2g'][even_chr_idx], 1, color='r',
        edgecolor='r', linewidth=0.1)

    # modify ticks
    ax3.set_xticks(tick_index)
    ax3.set_xticklabels(range(1,num_loci+1))
    ax3.tick_params(axis=u'both', which=u'both',length=0)

    # tuning figure
    ax3.set_xlim(0, num_loci)
    ax3.set_ylabel(trait2, fontsize=18)
    ax3.set_ylim(0.0, 0.002)
    ax3.yaxis.set_major_locator(ticker.MultipleLocator(0.001))
    plt.ticklabel_format(style='plain', axis='y', scilimits=(0,0))
    # make bar plot for h2g ##################################################
    

    # save figure ###########################################################
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%g'))
    plt.savefig(args.out, bbox_inches='tight')

# get command line
def get_command_line():
    
    parser = argparse.ArgumentParser(description='Make manhattan plot '\
        'for local genetic covariance estimates')
    
    parser.add_argument('--local-rhog-est', dest='local_rhog_file', type=str,
        help='Local genetic covariance estimation results', required=True)
    
    parser.add_argument('--local-hsqg-est', dest='local_hsqg_file', type=str,
        help='Local SNP-heritability estimation results',
        nargs=2, required=True)
    
    parser.add_argument('--trait-names', dest='trait_names', nargs=2, type=str,
        help='Names of traits', required=True)
    
    parser.add_argument('--out', type=str, dest='out', required=True,
        help='Output file name')

    args = parser.parse_args()
    
    return args

if(__name__ == '__main__'):
    main()
