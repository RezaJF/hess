#!/usr/bin/python
# (c) 2016-2021 Huwenbo Shi

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import argparse, math, sys
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
import matplotlib.gridspec as gridspec

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['ytick.labelsize'] = 16
matplotlib.rcParams['xtick.labelsize'] = 16
matplotlib.rc('font', family='sans-serif') 

def main():

    # get command line input
    args = get_command_line()

    # load local estimates file
    local_rhog = pd.read_table(args.local_rhog_file)
    local_hsqg1 = pd.read_table(args.local_hsqg_file[0])
    local_hsqg2 = pd.read_table(args.local_hsqg_file[1])

    # load gwas loci file
    gwas_loci1 = pd.read_table(args.gwas_loci_file[0], delim_whitespace=True)
    gwas_loci2 = pd.read_table(args.gwas_loci_file[1], delim_whitespace=True)

    # get trait names
    trait1, trait2 = args.traits
    trait1 = trait1.upper()
    trait2 = trait2.upper()

    # get loci
    loci = get_loci(local_rhog)
    mark_gwas1 = get_gwas_hit_loci(gwas_loci1, loci)
    mark_gwas2 = get_gwas_hit_loci(gwas_loci2, loci)
    mark_inter = get_intersection(mark_gwas1, mark_gwas2)
    mark_gwas1_spec = get_specific(mark_inter, mark_gwas1)
    mark_gwas2_spec = get_specific(mark_inter, mark_gwas2)
    mark_neither = get_neither(loci, mark_gwas1, mark_gwas2)

    # get gcor at ascertained loci
    gcor1, gcor_se1 = get_ascertained_gcor(local_rhog, local_hsqg1,
        local_hsqg2, mark_gwas1_spec)
    gcor2, gcor_se2 = get_ascertained_gcor(local_rhog, local_hsqg1,
        local_hsqg2, mark_gwas2_spec)
    gcor_int, gcor_se_int = get_ascertained_gcor(local_rhog, local_hsqg1,
        local_hsqg2, mark_inter)
    gcor_non, gcor_se_non = get_ascertained_gcor(local_rhog, local_hsqg1,
        local_hsqg2, mark_neither)
    
    # identify causality
    causal_trait = None
    effect_trait = None
    direction = None
    rng1 = (gcor1-gcor_se1*1.96, gcor1+gcor_se1*1.96)
    rng2 = (gcor2-gcor_se2*1.96, gcor2+gcor_se2*1.96)

    # overlaping
    if((rng1[1] > rng2[0] and rng1[1] < rng2[1]) or
       (rng2[1] > rng1[0] and rng2[1] < rng1[1])):
        pass
    # non-overlapping
    else:
        # only one of them overlapping 0?
        if((0 > rng1[0] and 0 < rng1[1]) and
           (not (0 > rng2[0] and 0 < rng2[1]))):
            causal_trait = trait2
            effect_trait = trait1
            if gcor2 > 0:
                direction = '\uparrow'
            else:
                direction = '\downarrow'
        elif((0 > rng2[0] and 0 < rng2[1]) and
           (not (0 > rng1[0] and 0 < rng1[1]))):
            causal_trait = trait1
            effect_trait = trait2
            if gcor1 > 0:
                direction = '\uparrow'
            else:
                direction = '\downarrow'

    # plot genetic correlation bar plot
    xvals = np.array([1,2,3,4])
    yvals = np.array([gcor1, gcor2, gcor_int, gcor_non])
    yerrs = np.array([gcor_se1, gcor_se2, gcor_se_int, gcor_se_non])
    outfile = open(args.out + '.log', 'w')
    outfile.write('trait1\tlocal_gcor1\tlocal_gcor1_se\t')
    outfile.write('trait2\tlocal_gcor2\tlocal_gcor2_se\n')
    outline = '%s\t%.2g\t%.2g\t%s\t%.2g\t%.2g\n' % (trait1, gcor1, gcor_se1,
        trait2, gcor2, gcor_se2)
    outfile.write(outline)
    outfile.close()
    fig, ax1 = plt.subplots(figsize=(4,6))
    fig.canvas.draw()
    wd = 0.5
    ax1.bar(xvals[0], yvals[0], width=wd, color='b',
        yerr=1.96*yerrs[0], ecolor='k', alpha=0.5, linewidth=0.5)
    ax1.bar(xvals[1], yvals[1], width=wd, color='r',
        yerr=1.96*yerrs[1], ecolor='k', alpha=0.5, linewidth=0.5)
    ax1.bar(xvals[2], yvals[2], width=wd, color='g',
        yerr=1.96*yerrs[2], ecolor='k', alpha=0.5, linewidth=0.5)
    ax1.bar(xvals[3], yvals[3], width=wd, color='g',
        yerr=1.96*yerrs[3], ecolor='k', alpha=0.5, linewidth=0.5)
    ax1.bar(xvals[0], yvals[0], width=wd, facecolor='None',
        yerr=1.96*yerrs[0], ecolor='k', alpha=1.0, linewidth=0.5)
    ax1.bar(xvals[1], yvals[1], width=wd, facecolor='None',
        yerr=1.96*yerrs[1], ecolor='k', alpha=1.0, linewidth=0.5)
    ax1.bar(xvals[2], yvals[2], width=wd, facecolor='None',
        yerr=1.96*yerrs[2], ecolor='k', alpha=1.0, linewidth=0.5)
    ax1.bar(xvals[3], yvals[3], width=wd, facecolor='None',
        yerr=1.96*yerrs[3], ecolor='k', alpha=1.0, linewidth=0.5)
    ax1.plot([0.5,1,2,3,4,4.5], [0, 0, 0, 0, 0, 0], "k--")
    ax1.set_xlim([0.75, 4.75])
    ax1.set_ylabel('local genetic correlation', fontsize=22)
    ax1.set_xticks([1,2,3,4])
    xlabels = [trait1+' specific', trait2+' specific',
               'intersection', 'neither']
    ax1.set_xticklabels(xlabels, rotation=45, ha='right', fontsize=16)

    title_str = ''
    if causal_trait != None and direction != None and effect_trait != None:
        title_str = r"%s $%s$ %s"% (causal_trait, direction, effect_trait)
    ax1.set_title(title_str, fontsize=22)
    
    # save the figure
    fig.subplots_adjust(left=0.35,bottom=0.3)
    plt.savefig(args.out)

# get loci
def get_loci(local_rhog):
    """
    Returns a map between chromosome number of a list of loci in the
    chromosome
    """
    loci = dict()
    for i in xrange(local_rhog.shape[0]):
        chrom = local_rhog['chr'][i]
        start, end = local_rhog['start'][i], local_rhog['end'][i]
        if(chrom not in loci): loci[chrom] = []
        loci[chrom].append((start, end))
    return loci

# get gwas hit loci
def get_gwas_hit_loci(gwas, loci):
    """
    Find the loci that harbors a GWAS hit
    """
    mark_gwas = set()

    # iterate through the chromosomes
    for chrom in xrange(1,23):
        loci_chrom = loci[chrom]

        # iterate through the loci in the chromosome
        for entry in loci_chrom:
            start = entry[0]
            end = entry[1]
            has_hit = False

            # get GWAS hit on the chromosome
            gwas_chrom = gwas[(gwas['CHR'] == chrom) & (gwas['BP'] >= start) &\
                              (gwas['BP'] <= end)]
            if(gwas_chrom.shape[0] != 0):
                mark_gwas.add('%d:%d-%d' % (chrom, start, end))
    
    return mark_gwas

# get intersection
def get_intersection(mark_gwas1, mark_gwas2):
    mark_inter = set()
    for loci in mark_gwas1:
        if(loci in mark_gwas2):
            mark_inter.add(loci)
    return mark_inter

# get specific
def get_specific(mark_inter, mark_gwas):
    spec = set()
    for loci in mark_gwas:
        if(loci not in mark_inter):
            spec.add(loci)
    return spec

# get neither
def get_neither(loci, mark_gwas1, mark_gwas2):
    neither = set()
    for chrom in xrange(1,23):
        loci_chrom = loci[chrom]
        for entry in loci_chrom:
            start = entry[0]
            end = entry[1]
            locus = '%d:%d-%d' % (chrom, start, end)
            if(locus not in mark_gwas1 and locus not in mark_gwas2):
                neither.add(locus)
    return neither

# get gcor at ascertained loci
def get_ascertained_gcor(gcov, hsq1, hsq2, mark_loci):
    
    # extract local rhog and local hsqg at the marked loci
    gcov_gwas_loci_list = []
    hsq1_gwas_loci_list = []
    hsq2_gwas_loci_list = []
    gcov_gwas_loci = 0.0
    for i in xrange(gcov.shape[0]):
        keyval = '{}:{}-{}'.format(gcov['chr'][i],
            gcov['start'][i], gcov['end'][i])
        if(keyval in mark_loci):
            gcov_gwas_loci += gcov['local_rhog'][i]
            hsq1_gwas_loci_list.append(hsq1['local_h2g'][i])
            hsq2_gwas_loci_list.append(hsq2['local_h2g'][i])
            gcov_gwas_loci_list.append(gcov['local_rhog'][i])
    gcov_gwas_loci_list = np.array(gcov_gwas_loci_list)
    hsq1_gwas_loci_list = np.array(hsq1_gwas_loci_list)
    hsq2_gwas_loci_list = np.array(hsq2_gwas_loci_list)
    
    # estimate gcor using all loci
    all_gcov = np.sum(gcov_gwas_loci_list)
    all_hsq1 = np.sum(hsq1_gwas_loci_list)
    all_hsq2 = np.sum(hsq2_gwas_loci_list)
    gcor_all_data = all_gcov/(np.sqrt(all_hsq1)*np.sqrt(all_hsq2))
    
    # get jackknife standard error
    gcor_var = 0.0
    for i in xrange(gcov_gwas_loci_list.shape[0]):
        tmp_gcov = all_gcov - gcov_gwas_loci_list[i]
        tmp_hsq1 = all_hsq1 - hsq1_gwas_loci_list[i]
        tmp_hsq2 = all_hsq2 - hsq2_gwas_loci_list[i]
        tmp_gcor = tmp_gcov/(np.sqrt(tmp_hsq1)*np.sqrt(tmp_hsq2))
        gcor_var += (tmp_gcor-gcor_all_data)**2.0
    gcor_se = np.sqrt(gcor_var)
    
    return (gcor_all_data, gcor_se)

# get command line
def get_command_line():

    parser = argparse.ArgumentParser(description='Infer putative causality' \
        ' between a pair of complex traits')
    
    parser.add_argument('--local-rhog-est', dest='local_rhog_file', type=str,
        help='Local genetic covariance estimates', required=True)
    
    parser.add_argument('--gwas-loci', dest='gwas_loci_file', type=str,
        help='GWAS hits files', required=True, nargs=2)
    
    parser.add_argument('--local-hsqg-est', dest='local_hsqg_file', type=str,
        help='Local SNP-heritability estimates', required=True, nargs=2)
    
    parser.add_argument('--trait-names', dest='traits', nargs=2, type=str,
        help='Names of the traits', required=True)
   
    parser.add_argument('--out', dest='out', type=str,
        help='Output file name', required=True)

    args = parser.parse_args()
    
    return args

if(__name__ == '__main__'):
    main()
