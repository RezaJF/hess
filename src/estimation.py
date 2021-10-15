import tarfile
import numpy as np
import pandas as pd
import scipy, scipy.stats
from sumstats import *
from refpanel import *

eps = 10.0**-8

def get_ld(snpdata):
    """
    Compute the LD matrix
    """
    ld = np.corrcoef(snpdata)
    ld = np.nan_to_num(ld)
    return ld

def eig_decomp(ld):
    """
    Perform eigenvalue decomposition on the LD matrix
    """
    ld_w, ld_v = np.linalg.eigh(ld)
    idx = ld_w.argsort()[::-1]
    ld_w = ld_w[idx]
    ld_v = ld_v[:,idx]
    return (ld_w.T, ld_v)

def local_hsqg_step1(refpanel_fnm, sumstats_fnm, partition_fnm,
    chrom, min_maf, out_fnm):
    """
    Implements the first step of local SNP-heritability estimation. Perform
    eigenvalue decomposition, compute the projection square.
    """

    # create files to write
    out_info = gzip.open('{}_chr{}.info.gz'.format(out_fnm, chrom), 'w')
    out_eig = gzip.open('{}_chr{}.eig.gz'.format(out_fnm, chrom), 'w')
    out_prjsq = gzip.open('{}_chr{}.prjsq.gz'.format(out_fnm, chrom), 'w')

    # load all files into memory
    partition = load_partition(partition_fnm, chrom)
    refpanel = PlinkReader(refpanel_fnm)
    sumstats = SumStats(sumstats_fnm, chrom)
    sumstats.filter_sumstats(refpanel.get_map())

    # iterate through loci
    for i in xrange(partition.shape[0]):
        start = partition['start'][i]
        stop = partition['stop'][i]
        
        # extract data in the locus and intersect
        snpmap_locus, snpdata_locus = refpanel.get_locus(start, stop, min_maf)
        sumstats_locus = sumstats.get_locus(start, stop)
        shared_snp = sumstats_locus.merge(snpmap_locus, on='SNP')[['SNP']]
        sumstats_locus = sumstats_locus[sumstats_locus['SNP']\
            .isin(shared_snp['SNP'])]

        # compute the ld matrix
        snpidx = np.where(snpmap_locus['SNP'].isin(shared_snp['SNP']))[0]
        snpdata_locus = snpdata_locus[snpidx, :]
        nsnp_locus = snpdata_locus.shape[0]
        logging.info('{} SNPs in locus chr{}:{}-{}'.format(nsnp_locus,
            chrom, start, stop))
        ld_locus = get_ld(snpdata_locus)
        
        # write out the information required for step 1
        local_hsqg_step1_helper(out_info, out_eig, out_prjsq, start, stop,
            ld_locus, sumstats_locus)

    # close all files to write
    out_info.close(); out_eig.close(); out_prjsq.close()

def local_hsqg_step1_helper(info, eig, prjsq, start, stop, ld, sumstats):
    """
    Helper function for step 1 of estimating local SNP-heritability
    """

    # perform eigenvalue decomposition
    nsnp = ld.shape[0]
    if(nsnp == 0):
        info.write('{}\t{}\t{}\t{}\t{:.1f}\n'.format(start, stop, 0, 0, 0.0))
        eig.write('\n'); prjsq.write('\n')
        return

    nindv = np.mean(sumstats['N'])
    rank = np.linalg.matrix_rank(ld)
    ld_w, ld_v = eig_decomp(ld)

    # write the information of the locus
    info.write('{}\t{}\t{}\t{}\t{:.1f}\n'.format(start,stop,nsnp,rank,nindv))

    # write the eigen value
    eig.write('\t'.join([str(ld_w[i]) for i in xrange(rank)])+'\n')

    # write the squared projection
    all_prjsq = []
    for i in xrange(rank):
        beta = sumstats['Z'] / np.sqrt(sumstats['N'])
        prj = np.dot(beta.T, ld_v[:, i])
        all_prjsq.append(prj**2)
    prjsq.write('\t'.join([str(all_prjsq[i]) for i in xrange(rank)])+'\n')

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
            logging.error('Missing step 1 results for chromosome {}'\
                .format(chrom))
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
        logging.error('Step 1 results contain different number of loci')
        sys.exit(1)

    # print out debugging info
    nloci = info.shape[0]
    logging.info('Loaded results for {} loci from step 1'.format(nloci))
    
    return (info, eig, prjsq)

def local_hsqg_step2(prefix, max_num_eig, min_eigval, reinflate,
    gwse_thres, tot_hsq, out_fnm):
    """
    Implements the second step of HESS -- uses results from step 1, corrects
    for biases, and obtains standard error
    """

    # make sure reinflate is not a list
    if reinflate == None: reinflate = [1.0]
    if len(reinflate) > 1:
        logging.error('Specified lambda gc should be a single value')
        sys.exit(1)
    reinflate = reinflate[0]

    # load results from step 1
    info, eig, prjsq = load_local_hsqg_step1(prefix)

    # estimate local SNP-heritability without genome-wide SNP-heritability
    order = ['chr','start','end','num_snp','k','local_h2g','var','se','z','p']
    if tot_hsq == None:
        local_hsqg_est = local_hsqg_step2_helper(info, eig, prjsq, max_num_eig,
            min_eigval, reinflate, gwse_thres)
        local_hsqg_est = local_hsqg_est[order] 
        local_hsqg_est.to_csv(out_fnm+'.txt', '\t', header=True,
            index=None, float_format='%.5g')
    else:
        local_hsqg_est = local_hsqg_step2_helper_tot_hsq(info, eig, prjsq,
            max_num_eig, min_eigval, reinflate, gwse_thres, tot_hsq)
        local_hsqg_est = local_hsqg_est[order] 
        local_hsqg_est.to_csv(out_fnm+'.txt', '\t', header=True,
            index=None, float_format='%.5g')

def local_quad_form(eig, prjsq, max_k, min_eigval, reinflate):
    """
    Compute the quadratic form (beta_gwas' * LD_inv * beta_gwas) at each locus
    """
    
    # log the lambda gc to reinflate
    logging.info('Re-inflate the summary statistics with '\
        'lambda_gc: {:.4g}'.format(reinflate))

    all_sum = []; all_k = []
    nloci = len(eig)
    for i in xrange(nloci):
        k = min(max_k, np.where(eig[i] > min_eigval)[0].size)
        tmp = np.divide(prjsq[i][0:k], eig[i][0:k]+eps) * reinflate
        all_sum.append(np.sum(tmp))
        all_k.append(float(k))
    return pd.DataFrame({'sum': all_sum, 'k':all_k})

def local_hsqg_var(info, local_hsqg_est):
    """
    Estimate the variance of local SNP-heritability
    """

    # get total SNP-heritability estimate
    tot_hsqg = np.sum(local_hsqg_est)
    nloci = local_hsqg_est.shape[0]

    # estimate variance
    A = np.zeros((nloci, nloci))
    for i in xrange(nloci):
        A[i,:] = np.square(info['rank'] / (info['N'] - info['rank'] + eps))
    A[range(nloci), range(nloci)] = 1.0
    b = np.square(info['N'] / (info['N'] - info['rank'] + eps))
    b *= (2.0*info['rank']*((1.0-tot_hsqg)/(info['N']+eps))+4.0*local_hsqg_est)
    b *= ((1.0-tot_hsqg) / (info['N'] + eps))

    # check the rank of A
    rank = np.linalg.matrix_rank(A)
    if rank < nloci:
        logging.error('Rank of A less than the number of loci when estimating'\
            ' the variance of the local SNP-heritability estimates.')
        sys.exit(1)

    # get the variance estimates
    local_hsqg_est_var = np.dot(np.linalg.pinv(A), b)
    tot_hsqg_var = np.sum(local_hsqg_est_var)

    # check if the total variance is positive
    tot_hsqg_se = np.nan
    if tot_hsqg_var > 0.0:
        tot_hsqg_se = np.sqrt(tot_hsqg_var)
    else:
        logging.warning('Variance of total SNP-heritability estimates is'\
            ' negative. This is likely due to small sample size of the GWAS.')

    # log the total SNP-heritability and standard error
    logging.info('Total SNP-heritability estimate: {:.3g} ({:.3g})'\
        .format(tot_hsqg, tot_hsqg_se))

    return pd.Series(local_hsqg_est_var)

def local_hsqg_step2_helper(info, eig, prjsq, max_num_eig, min_eigval,
    reinflate, gwse_thres):
    """
    Estimate local SNP-heritability when total SNP-heritability is not
    provided
    """

    # find the max number of eigenvectors to use based on the sample size
    if gwse_thres == None: gwse_thres = 2.0
    nloci = info.shape[0]; nsnp = np.sum(info['nsnp'])
    nindv = np.sum(info['N']*info['nsnp']) / nsnp
    max_k = (gwse_thres * nindv - nindv) / (gwse_thres * nloci + eps)
    max_k = min(int(np.ceil(max_k)), max_num_eig)
    
    # log the number of SNPs and sample size
    logging.info('Using {} SNPs with average sample size {}'\
        .format(nsnp, nindv))

    # get the quadratic form
    quad_form = local_quad_form(eig, prjsq, max_k, min_eigval, reinflate)
    
    # correct for bias in the quadratic form
    A = np.diag(info['N'])
    for i in xrange(nloci):
        A[i,:] -= quad_form['k']
    b = info['N'] * quad_form['sum'] - quad_form['k']
    
    # check the rank of A
    rank = np.linalg.matrix_rank(A)
    if rank < nloci:
        logging.error('Rank of A less than the number of loci. There' \
            ' might be loci with no SNP.')
        sys.exit(1)

    # obtain the bias-corrected estimates and variance
    local_hsqg_est = np.dot(np.linalg.pinv(A), b)
    local_hsqg_est_var = local_hsqg_var(info, local_hsqg_est)
    
    # obtain test statistics
    se = np.sqrt(local_hsqg_est_var.astype(np.float32))
    zsc = local_hsqg_est / se
    zsc[(zsc < 0) | zsc.isnull()] = 0.0
    pval = scipy.stats.norm.sf(np.fabs(zsc))

    # construct the result data frame
    result = pd.DataFrame({'chr': info['CHR'], 'start': info['start'],
        'end': info['stop'], 'num_snp': info['nsnp'],
        'k': quad_form['k'].astype(np.int), 'local_h2g':local_hsqg_est,
        'var': local_hsqg_est_var, 'se': se, 'z': zsc, 'p': pval})

    return result

def local_hsqg_step2_helper_tot_hsq(info, eig, prjsq, max_num_eig, min_eigval,
    reinflate, gwse_thres, tot_hsq):
    """
    Estimate local SNP-heritability when total SNP-heritability and its
    standard eror are provided
    """
    
    # find the max number of eigenvectors to use based on the sample size
    if gwse_thres == None: gwse_thres = 0.5
    nloci = info.shape[0]; nsnp = np.sum(info['nsnp'])
    nindv = np.sum(info['N']*info['nsnp']) / nsnp
    max_k = (gwse_thres * nindv) / (nloci + eps)
    max_k = min(int(np.ceil(max_k)), max_num_eig)
    
    # log the number of SNPs and sample size
    logging.info('Using {} SNPs with average sample size {}'\
        .format(nsnp, nindv))

    # get the quadratic form
    tot_hsq_val, tot_hsqg_se = tot_hsq
    quad_form = local_quad_form(eig, prjsq, max_k, min_eigval, reinflate)

    # correct for bias in the quadratic form
    bias = (1.0-tot_hsq_val)*quad_form['k']/(info['N']+eps)
    local_hsqg_est = quad_form['sum'] - bias
    tot_hsqg = np.sum(local_hsqg_est)

    # get variance estimate
    local_hsqg_est_var = np.square(info['N'] / (info['N']-info['rank']+eps))
    local_hsqg_est_var *= (2.0*info['rank']*((1.0-local_hsqg_est)\
        /(info['N']+eps))+4.0*local_hsqg_est)
    local_hsqg_est_var *= ((1.0 - local_hsqg_est) / (info['N'] + eps))
    local_hsqg_est_var += np.square(tot_hsqg_se*info['rank']/(info['N']+eps))
    tot_hsqg_var = np.sum(local_hsqg_est_var)

    # check if the total variance is positive
    tot_hsqg_se = np.nan
    if tot_hsqg_var > 0.0:
        tot_hsqg_se = np.sqrt(tot_hsqg_var)
    else:
        logging.warning('Variance of total SNP-heritability estimates is'\
            ' negative. This is likely due to small sample size of the GWAS.')

    # log the total SNP-heritability and standard error
    logging.info('Total SNP-heritability estimate: {:.3g} ({:.3g})'\
        .format(tot_hsqg, tot_hsqg_se))
  
    # obtain test statistics
    se = np.sqrt(local_hsqg_est_var)
    zsc = local_hsqg_est / se
    zsc[(zsc < 0) | zsc.isnull()] = 0.0
    pval = scipy.stats.norm.sf(np.fabs(zsc))

    # construct the result data frame
    result = pd.DataFrame({'chr': info['CHR'], 'start': info['start'],
        'end': info['stop'], 'num_snp': info['nsnp'],
        'k': quad_form['k'].astype(np.int), 'local_h2g':local_hsqg_est,
        'var': local_hsqg_est_var, 'se': se, 'z': zsc, 'p': pval})

    return result

###############################################################################

def local_rhog_step1(refpanel_fnm, sumstats_fnm, partition_fnm,
    chrom, min_maf, out_fnm):
    """
    Implements the first step of local genetic covariance estimation. Perform
    eigenvalue decomposition, compute the projection square of each trait,
    and projection product across the traits.
    """

    # make sure that two summary stats are given
    if type(sumstats_fnm) != list or len(sumstats_fnm) != 2:
        logging.error('Must provide two summary statistics data')
        sys.exit(1)

    # create files to write
    out_info1 = gzip.open('{}_trait1_chr{}.info.gz'.format(out_fnm, chrom),'w')
    out_eig1 = gzip.open('{}_trait1_chr{}.eig.gz'.format(out_fnm, chrom), 'w')
    out_prjsq1=gzip.open('{}_trait1_chr{}.prjsq.gz'.format(out_fnm,chrom),'w')

    out_info2 = gzip.open('{}_trait2_chr{}.info.gz'.format(out_fnm, chrom),'w')
    out_eig2 = gzip.open('{}_trait2_chr{}.eig.gz'.format(out_fnm, chrom), 'w')
    out_prjsq2=gzip.open('{}_trait2_chr{}.prjsq.gz'.format(out_fnm,chrom),'w')

    out_eig = gzip.open('{}_chr{}.eig.gz'.format(out_fnm, chrom), 'w')
    out_prjprod = gzip.open('{}_chr{}.prjprod.gz'.format(out_fnm, chrom), 'w')

    # load all files into memory
    partition = load_partition(partition_fnm, chrom)
    refpanel = PlinkReader(refpanel_fnm)
    sumstats1 = SumStats(sumstats_fnm[0], chrom)
    sumstats1.filter_sumstats(refpanel.get_map())
    sumstats2 = SumStats(sumstats_fnm[1], chrom)
    sumstats2.filter_sumstats(refpanel.get_map())

    # iterate through loci
    for i in xrange(partition.shape[0]):
        start = partition['start'][i]
        stop = partition['stop'][i]
       
        # extract data in the locus and intersect
        snpmap_locus, snpdata_locus = refpanel.get_locus(start, stop, min_maf)
        sumstats1_locus = sumstats1.get_locus(start, stop)
        sumstats2_locus = sumstats2.get_locus(start, stop)
        shared_snp = sumstats1_locus.merge(snpmap_locus, on='SNP')[['SNP']]
        shared_snp = sumstats2_locus.merge(shared_snp, on='SNP')[['SNP']]
        sumstats1_locus = sumstats1_locus[sumstats1_locus['SNP']\
            .isin(shared_snp['SNP'])]
        sumstats2_locus = sumstats2_locus[sumstats2_locus['SNP']\
            .isin(shared_snp['SNP'])]

        # compute the ld matrix
        snpidx = np.where(snpmap_locus['SNP'].isin(shared_snp['SNP']))[0]
        snpdata_locus = snpdata_locus[snpidx, :]
        nsnp_locus = snpdata_locus.shape[0]
        logging.info('{} SNPs in locus chr{}:{}-{}'.format(nsnp_locus,
            chrom, start, stop))
        ld_locus = get_ld(snpdata_locus)
        
        # write out the information required for step 1
        local_hsqg_step1_helper(out_info1, out_eig1, out_prjsq1, start, stop,
            ld_locus, sumstats1_locus)
        local_hsqg_step1_helper(out_info2, out_eig2, out_prjsq2, start, stop,
            ld_locus, sumstats2_locus)
        local_rhog_step1_helper(out_eig, out_prjprod, ld_locus,
            sumstats1_locus, sumstats2_locus)

    # close all files to write
    out_info1.close(); out_eig1.close(); out_prjsq1.close()
    out_info2.close(); out_eig2.close(); out_prjsq2.close()
    out_eig.close(); out_prjprod.close()

def local_rhog_step1_helper(eig, prjprod, ld, sumstats1, sumstats2):
    """
    Helper function for step 1 of estimating local SNP-heritability
    """

    # perform eigenvalue decomposition
    nsnp = ld.shape[0]
    if(nsnp == 0):
        eig.write('\n'); prjprod.write('\n')
        return
    
    rank = np.linalg.matrix_rank(ld)
    ld_w, ld_v = eig_decomp(ld)

    # write the eigen value
    eig.write('\t'.join([str(ld_w[i]) for i in xrange(rank)])+'\n')

    # write the squared projection
    all_prjprod = []
    for i in xrange(rank):
        beta = sumstats1['Z'] / np.sqrt(sumstats1['N'])
        gamma = sumstats2['Z'] / np.sqrt(sumstats2['N'])
        prj1 = np.dot(beta.T, ld_v[:, i])
        prj2 = np.dot(gamma.T, ld_v[:, i])
        all_prjprod.append(prj1*prj2)
    prjprod.write('\t'.join([str(all_prjprod[i]) for i in xrange(rank)])+'\n')

def local_rhog_step2(prefix, local_hsqg_est_fnm, max_num_eig, min_eigval,
    reinflate, gwse_thres, pheno_cor, num_shared, out_fnm):
    """
    Implements step 2 for estimating local genetic covariance
    """

    # make sure reinflate is a list
    if reinflate == None: reinflate = [1.0, 1.0]
    if type(reinflate) != list or len(reinflate) != 2:
        logging.error('Must provide two lambda gc')
        sys.exit(1)

    # load the local SNP-heritability estimates
    local_hsqg_est1, local_hsqg_est2 = load_local_hsqg_est(local_hsqg_est_fnm)

    # load result from step 1
    info1, info2, eig, prjprod = load_local_rhog_step1(prefix)

    # compute the bilinear form
    if gwse_thres == None: gwse_thres = 2.0
    nloci = info1.shape[0]; nsnp = np.sum(info1['nsnp'])
    nindv1 = np.sum(info1['N']*info1['nsnp']) / nsnp
    nindv2 = np.sum(info2['N']*info2['nsnp']) / nsnp
    max_k = max_num_eig
    if(num_shared != 0):
        max_k = (gwse_thres * nindv1 * nindv2 - nindv1 * nindv2)
        max_k /= (gwse_thres * num_shared * nloci + eps)
    max_k = min(int(np.ceil(max_k)), max_num_eig)
    
    logging.info('Using {} SNPs with average sample size {} for trait 1'\
        .format(nsnp, nindv1))
    logging.info('Using {} SNPs with average sample size {} for trait 2'\
        .format(nsnp, nindv2))
    
    bilin_form = local_bilin_form(eig, prjprod, max_k, min_eigval, reinflate)

    logging.info('Using phenotypic correlation: {:.3g}'.format(pheno_cor))

    # correct for bias in the bilinar form
    A = np.diag(info1['N']*info2['N'])
    for i in xrange(nloci):
        A[i,:] -= num_shared * bilin_form['k']
    b = info1['N']*info2['N']*bilin_form['sum'] - \
        num_shared*pheno_cor*bilin_form['k']
    
    # check the rank of A
    rank = np.linalg.matrix_rank(A)
    if rank < nloci:
        logging.error('Rank of A less than the number of loci. There' \
            ' might be loci with no SNP.')
        sys.exit(1)

    # obtain the bias-corrected estimates
    local_rhog_est = np.dot(np.linalg.pinv(A), b)
    tot_rhog = np.sum(local_rhog_est)

    # obtain the standard error
    tot_hsqg1 = np.sum(local_hsqg_est1['local_h2g'])
    tot_hsqg2 = np.sum(local_hsqg_est2['local_h2g'])
    rhoe = pheno_cor-tot_rhog
    if num_shared == 0: rhoe = 0.0
    sige1 = 1.0-tot_hsqg1; sige2 = 1.0-tot_hsqg2
    
    A = np.zeros((nloci, nloci))
    n1n2 = info1['N']*info2['N']
    nsp = num_shared*info1['rank']
    for i in xrange(nloci):
        A[i,:] = -1.0 * np.square(nsp / (n1n2 - nsp + eps))
    A[range(nloci), range(nloci)] = 1.0
    b = np.square(n1n2 / (n1n2 -nsp + eps))
    term1 = np.square(rhoe*info1['rank']*info1['rank'] / (n1n2 + eps))
    term2 = (sige1**2.0)*(sige2**2.0)*info1['rank'] / (n1n2 + eps)
    term3 = sige2*local_hsqg_est1['local_h2g'] / (info2['N'] + eps)
    term4 = sige1*local_hsqg_est2['local_h2g'] / (info1['N'] + eps)
    term5 = 2.0*num_shared*local_rhog_est*tot_rhog / (n1n2 + eps)
    b = b * (term1 + term2 + term3 + term4 + term5)

    # make sure A is invertible
    rank = np.linalg.matrix_rank(A)
    if rank < nloci:
        logging.error('Rank of A less than the number of loci when' \
            ' quantifying the variance estimates')
        sys.exit(1)

    local_rhog_est_var = np.dot(np.linalg.pinv(A), b)
    local_rhog_est_var = pd.Series(local_rhog_est_var)
    tot_rhog_var = np.sum(local_rhog_est_var)

    # check if the total variance is positive
    tot_rhog_se = np.nan
    if tot_rhog_var > 0.0:
        tot_rhog_se = np.sqrt(tot_rhog_var)
    else:
        logging.warning('Variance of genome-wide genetic covariance estimate'\
            ' is negative. This is likely due to small sample size of'\
            ' the GWAS.')

    # log the genome-wide genetic covariance and standard error
    logging.info('Genome-wide genetic covariance estimate: {:.3g} ({:.3g})'\
        .format(tot_rhog, tot_rhog_se))

    # obtain test statistics
    se = np.sqrt(local_rhog_est_var)
    zsc = local_rhog_est / se
    zsc[zsc.isnull()] = 0.0
    pval = scipy.stats.norm.sf(np.fabs(zsc))*2.0

    # construct the result data frame
    result = pd.DataFrame({'chr': info1['CHR'], 'start': info1['start'],
        'end': info1['stop'], 'num_snp': info1['nsnp'],
        'k': bilin_form['k'].astype(np.int), 'local_rhog':local_rhog_est,
        'var': local_rhog_est_var, 'se':se, 'z': zsc, 'p': pval})

    result = result[['chr','start','end','num_snp','k','local_rhog',
        'var', 'se', 'z', 'p']]
    result.to_csv(out_fnm+'.txt', '\t', header=True,
        index=None, float_format='%.5g')

    # estimate genome-wide genetic correlation
    estimate_gw_rg(result, local_hsqg_est1, local_hsqg_est2)

def estimate_gw_rg(local_rhog_est, local_hsqg_est1, local_hsqg_est2):
    """
    Estimate genome-wide genetic correlation and its standard error
    """
    
    # estimate using all data
    tot_rhog = np.sum(local_rhog_est['local_rhog'])
    tot_hsqg1 = np.sum(local_hsqg_est1['local_h2g'])
    tot_hsqg2 = np.sum(local_hsqg_est2['local_h2g'])

    if tot_hsqg1 < 0 or tot_hsqg2 < 0:
        logging.warning('Negative total SNP-heritability')
        gcor_est = 0.0
    else:
        gcor_est = tot_rhog / np.sqrt(tot_hsqg1 * tot_hsqg2)

    # get jackknife se
    jk_est = []; nloci = local_rhog_est.shape[0]
    for i in xrange(nloci):
        local_rhog_est_jk = local_rhog_est.drop(i)
        local_hsqg_est_jk1 = local_hsqg_est1.drop(i)
        local_hsqg_est_jk2 = local_hsqg_est2.drop(i)
        
        tot_rhog_jk = np.sum(local_rhog_est_jk['local_rhog'])
        tot_hsqg_jk1 = np.sum(local_hsqg_est_jk1['local_h2g'])
        tot_hsqg_jk2 = np.sum(local_hsqg_est_jk2['local_h2g'])

        if tot_hsqg_jk1 < 0 or tot_hsqg_jk2 < 0: jk_est.append(0.0)
        else: jk_est.append(tot_rhog_jk/np.sqrt(tot_hsqg_jk1*tot_hsqg_jk2))

    jk_est = np.array(jk_est)
    jk_se = np.sqrt((nloci-1.0)*np.mean(np.square(jk_est-gcor_est)))

    logging.info('Genome-wide genetic correlation estimate: {:.3g} ({:.3g})'\
        .format(gcor_est, jk_se))

def load_local_hsqg_est(local_hsqg_est_fnm):
    """
    Load the local SNP-heritability estimates
    """

    # sanity check of the input
    if type(local_hsqg_est_fnm) != list or len(local_hsqg_est_fnm) != 2:
        logging.error('Must provide two local SNP-heritability estimates')
        sys.exit(1)

    # load local SNP-heritability estimates
    local_hsqg_est1=pd.read_table(local_hsqg_est_fnm[0],delim_whitespace=True)
    local_hsqg_est2=pd.read_table(local_hsqg_est_fnm[1],delim_whitespace=True)
    
    # make sure they have the same number of loci
    if local_hsqg_est1.shape[0] != local_hsqg_est2.shape[0]:
        logging.error('Local SNP-heritability estimates have different number'
            'of loci')
        sys.exit(1)
    
    # log the total SNP-heritability
    tot_hsqg1 = np.sum(local_hsqg_est1['local_h2g'])
    tot_hsqg2 = np.sum(local_hsqg_est2['local_h2g'])
    tot_hsqg1_var = np.sum(local_hsqg_est1['var'])
    tot_hsqg2_var = np.sum(local_hsqg_est2['var'])

    # give a warning if variance is less than 0
    if tot_hsqg1_var < 0 or tot_hsqg2_var < 0:
        logging.warning('Negative variance of total SNP-heritability detected')

    tot_hsqg1_se = np.nan; tot_hsqg2_se = np.nan
    if tot_hsqg1_var > 0: tot_hsqg1_se = np.sqrt(tot_hsqg1_var)
    if tot_hsqg2_var > 0: tot_hsqg2_se = np.sqrt(tot_hsqg2_var)
    logging.info('Total SNP-heritability of trait 1: {:.3g} ({:.3g})'\
        .format(tot_hsqg1, tot_hsqg1_se))
    logging.info('Total SNP-heritability of trait 2: {:.3g} ({:.3g})'\
        .format(tot_hsqg2, tot_hsqg2_se))

    return (local_hsqg_est1, local_hsqg_est2)

def local_bilin_form(eig, prjprod, max_k, min_eigval, reinflate):
    """
    Compute the bilinear form (beta_gwas' * LD_inv * gamma_gwas) at each locus
    """
    
    # log the lambda gc to reinflate
    logging.info('Re-inflate the summary statistics with '\
        'lambda_gc: {:.4g} and {:.4g}'.format(reinflate[0], reinflate[1]))

    all_sum = []; all_k = []
    nloci = len(eig)
    for i in xrange(nloci):
        k = min(max_k, np.where(eig[i] > min_eigval)[0].size)
        tmp = np.divide(prjprod[i][0:k], eig[i][0:k]+eps)
        tmp *= np.sqrt(reinflate[0] * reinflate[1])
        all_sum.append(np.sum(tmp))
        all_k.append(float(k))
    return pd.DataFrame({'sum': all_sum, 'k':all_k})


def load_local_rhog_step1(prefix):
    """
    Load results from step 1 for estimating local genetic covariance
    """

    # check if all files exist
    for chrom in xrange(1, 23):
        info1_f = '{}_trait1_chr{}.info.gz'.format(prefix, chrom)
        info2_f = '{}_trait2_chr{}.info.gz'.format(prefix, chrom)
        eig_f = '{}_chr{}.eig.gz'.format(prefix, chrom)
        prjprod_f = '{}_chr{}.prjprod.gz'.format(prefix, chrom)
        if((not os.path.exists(info1_f)) or (not os.path.exists(info2_f)) or
           (not os.path.exists(eig_f)) or (not os.path.exists(prjprod_f))):
            logging.error('Missing step 1 results for chromosome {}'\
                .format(chrom))
            sys.exit(1)

    # iterate through chromosomes
    eig = []; prjprod = []
    info1 = pd.DataFrame(); info2 = pd.DataFrame()
    for chrom in xrange(1, 23):
        
        # load the information about each locus into memory
        info1_f = '{}_trait1_chr{}.info.gz'.format(prefix, chrom)
        info1_chr = pd.read_table(info1_f, delim_whitespace=True, header=None,
            compression='gzip', names=['start', 'stop', 'nsnp', 'rank', 'N'])
        info1_chr['CHR'] = chrom
        info1 = pd.concat((info1, info1_chr), axis=0)

        info2_f = '{}_trait2_chr{}.info.gz'.format(prefix, chrom)
        info2_chr = pd.read_table(info2_f, delim_whitespace=True, header=None,
            compression='gzip', names=['start', 'stop', 'nsnp', 'rank', 'N'])
        info2_chr['CHR'] = chrom
        info2 = pd.concat((info2, info2_chr), axis=0)

        # load the eigenvalues and squared projections into memory
        eig_f = gzip.open('{}_chr{}.eig.gz'.format(prefix, chrom), 'r')
        prjprod_f = gzip.open('{}_chr{}.prjprod.gz'.format(prefix, chrom), 'r')
        for line in eig_f:
            eig.append(np.array(line.strip().split()).astype(np.float))
        for line in prjprod_f:
            prjprod.append(np.array(line.strip().split()).astype(np.float))
        eig_f.close(); prjprod_f.close()

    # reset the index of info
    info1 = info1.reset_index(drop=True)
    info2 = info2.reset_index(drop=True)

    # check if info, eig, and prjsq have the same length
    if(info1.shape != info2.shape or info2.shape[0] != len(eig) or
       len(eig) != len(prjprod)):
        logging.error('Step 1 results contain different number of loci')
        sys.exit(1)

    # print out debugging info
    nloci = info1.shape[0]
    logging.info('Loaded results for {} loci from step 1'.format(nloci))
    
    return (info1, info2, eig, prjprod)

