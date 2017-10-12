import numpy as np, numpy.linalg
import pandas as pd
import sys,argparse,os,gzip

eps = 10.0**-8

def main():
   
    # get command line
    args = get_command_line()

    # load step 1 results
    info, eig, prjsq = load_local_hsqg_step1(args.prefix)

    # compute emprical and theoretical
    nloci = info.shape[0]
    empirical = local_quad_form(info,eig,prjsq, args.num_eig, args.min_eigval)   
    empirical = empirical.sort_values(by=['emp']).reset_index(drop=True)
    empirical = empirical[0:int(nloci*args.pct)]
  
    # esimate the lambda gc
    emp = np.reshape(empirical['emp'].values, (-1,1))
    theo = empirical['theo'].values
    lambda_gc = max(np.linalg.lstsq(emp, theo)[0][0], 1.0)
    print lambda_gc

def load_local_hsqg_step1(prefix):
    """
    Load results from step 1 for estimating local SNP-heritability
    """

    # check if all files exist
    for chrom in xrange(1, 23):
        info_f = '{}_chr{}.info.gz'.format(prefix, chrom)
        eig_f = '{}_chr{}.eig.gz'.format(prefix, chrom)
        prjsq_f = '{}_chr{}.prjsq.gz'.format(prefix, chrom)
        if((not os.path.exists(info_f)) or (not os.path.exists(eig_f)) or \
           (not os.path.exists(prjsq_f))):
            sys.exit(1)

    # iterate through chromosomes
    info = pd.DataFrame(); eig = []; prjsq = []
    for chrom in xrange(1, 23):
        
        # load the information about each locus into memory
        info_f = '{}_chr{}.info.gz'.format(prefix, chrom)
        info_chr = pd.read_table(info_f, delim_whitespace=True, header=None,
            compression='gzip', names=['start', 'stop', 'nsnp', 'rank', 'N'])
        info_chr['CHR'] = chrom
        info = pd.concat((info, info_chr), axis=0)

        # load the eigenvalues and squared projections into memory
        eig_f = gzip.open('{}_chr{}.eig.gz'.format(prefix, chrom), 'r')
        prjsq_f = gzip.open('{}_chr{}.prjsq.gz'.format(prefix, chrom), 'r')
        for line in eig_f:
            eig.append(np.array(line.strip().split()).astype(np.float))
        for line in prjsq_f:
            prjsq.append(np.array(line.strip().split()).astype(np.float))
        eig_f.close(); prjsq_f.close()

    # reset the index of info
    info = info.reset_index(drop=True)

    # check if info, eig, and prjsq have the same length
    if info.shape[0] != len(eig) or len(eig) != len(prjsq):
        sys.exit(1)

    # print out debugging info
    nloci = info.shape[0]
    
    return (info, eig, prjsq)

def local_quad_form(info, eig, prjsq, max_k, min_eigval):
    """
    Compute the quadratic form (beta_gwas' * LD_inv * beta_gwas) at each locus
    """
    
    all_sum = []; all_k = []; all_theo = []
    nloci = len(eig)
    for i in xrange(nloci):
        k = min(max_k, np.where(eig[i] > min_eigval)[0].size)
        tmp = np.divide(prjsq[i][0:k], eig[i][0:k]+eps)
        all_sum.append(np.sum(tmp))
        all_k.append(float(k))
        all_theo.append(float(k)/info['N'][i])

    return pd.DataFrame({'emp': all_sum, 'k':all_k, 'theo': all_theo})

# get command line
def get_command_line():
    
    parser = argparse.ArgumentParser(description='Estimate lambda gc '
        'to re-inflate')

    parser.add_argument('--prefix', dest='prefix', type=str, required=True,
        help='Local SNP-heritability estimation step 1 output prefix')
   
    parser.add_argument('--num-eig', dest='num_eig', type=int, required=False,
        help='Number of eigenvectors to use', default=50)
   
    parser.add_argument('--pct', dest='pct', type=float, required=False,
        help='Percent of loci to use to estimate lambda gc', default=0.5)

    parser.add_argument('--min-eigval', dest='min_eigval', type=float,
        help='Minimum eigenvalue', default=1.0, required=False)
 
    args = parser.parse_args()
    
    return args

if(__name__ == '__main__'):
    main()
