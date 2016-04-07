#!/usr/bin/python

from optparse import OptionParser
import numpy as np, numpy.linalg
import argparse, math, sys

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
    snp_idx,all_snp_leg = load_legend(leg_file)

    # load zscore file
    snp_beta,all_snp_beta = load_beta(zsc_file)
    
    # load partition
    part = load_partition(part_file)
    
    # output eigenvalue and projection squared
    output_eig_prjsq(chrom, snp_idx, all_snp_leg, snp_beta,
        all_snp_beta, part, ref_file, out_file)

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

# standardize snps
def standardize_rows(mat):
    nrow = mat.shape[0]
    for i in xrange(nrow):
        mu = np.mean(mat[i,:])
        var = np.var(mat[i,:])
        std = math.sqrt(var)
        mat[i,:] = mat[i,:]-mu
        mat[i,:] = mat[i,:]/std
    return mat

# get eigenvalue decomposition
def eig_decomp(locus_ld):
    ld_w,ld_v = np.linalg.eigh(locus_ld)
    idx = ld_w.argsort()[::-1]
    ld_w = ld_w[idx]
    ld_v = ld_v[:,idx]
    ld_w = np.matrix(ld_w).transpose()
    ld_v = np.matrix(ld_v)
    return (ld_w, ld_v)

# get lines of the reference file to load
def get_load_line_idx(snp_idx, all_snp_leg):
    load_line_idx = set()
    for snp in all_snp_leg:
        if(snp[0] in snp_idx):
            load_line_idx.add(snp_idx[snp[0]])
    return load_line_idx

# get ld matrix
def get_ld(mat):
    ld = np.corrcoef(mat)
    ld = np.nan_to_num(ld)
    return np.matrix(ld)

# load haplotypes
def load_ref_file(ref_file, load_line_idx, all_snp_leg, line_idx):
    snp_data = dict()
    new_line_idx = line_idx
    num_load = 0
    while(True):
        line = ref_file.readline()
        if(not line):
            break
        if(new_line_idx in load_line_idx):
            snp = all_snp_leg[new_line_idx]
            line = line.strip()
            cols = line.split()
            snp_data[snp] = [float(cols[i]) for i in xrange(len(cols))]
            num_load += 1
        new_line_idx += 1
        if(num_load == len(load_line_idx)):
            break
    last_fp = ref_file.tell()
    return snp_data,new_line_idx,last_fp

# get haplotypes in a locus
def get_locus_refpanel(snp_hap, locus_snp):
    haps = []
    for i in xrange(len(locus_snp)):
        haps.append(snp_hap[locus_snp[i][0]])
    haps = np.matrix(haps)
    return haps[:,0::2]+haps[:,1::2]

# get beta for window snps
def get_locus_snp_beta(snp_beta, locus_snp):
    beta = []
    for snp in locus_snp:
        beta.append(snp_beta[snp[0]])
    return beta

# load legend
def load_legend(leg_file):
    all_snps = []
    leg_file = open(leg_file, 'r')
    flr = False
    snp_idx = dict()
    idx = 0
    for line in leg_file:
        line = line.strip()
        cols = line.split()
        snpid = cols[0]
        all_snps.append(snpid)
        snp_idx[snpid] = idx
        idx += 1
    leg_file.close()
    return (snp_idx, all_snps)

# load beta
def load_beta(zsc_file):
    all_snps = []
    snp_beta = dict()
    zsc_file = open(zsc_file, 'r')
    for line in zsc_file:
        line = line.strip()
        cols = line.split()
        beta = float(cols[4])/math.sqrt(float(cols[5]))
        snp_beta[cols[0]] = beta
        all_snps.append((cols[0], int(cols[1]), float(cols[5])))
    zsc_file.close()
    return (snp_beta, all_snps)

# make windows
def load_partition(part_file):
    part = []
    flr = False
    part_file = open(part_file, 'r')
    for line in part_file:
        if(flr == False):
            flr = True
            continue
        line = line.strip()
        cols = line.split()
        st = int(cols[1])
        ed = int(cols[2])-1
        part.append((st,ed))
    part_file.close()
    return part

# get legend of snps in the locus
# start_idx: index where the search in the legend should start
# return: legend for snps in the locus, and where the search ends
def get_locus_snp(all_snp, locus, start_idx):
    locus_snp = []
    end_idx = start_idx
    for i in xrange(start_idx, len(all_snp)):
        snp = all_snp[i]
        pos = int(snp[1])
        if(pos >= locus[0] and pos <= locus[1]):
            locus_snp.append(snp)
        if(pos > locus[1]):
            break
        end_idx += 1
    return (locus_snp,end_idx)

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
