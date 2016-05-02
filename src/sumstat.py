# (c) 2016-2021 Huwenbo Shi


import sys


"""
description:
    filter out any snp in z-score file that is not in reference panel
    filter out any snp with strand ambiguilty
    fix the sign of beta to be consistent with reference panel
arguments:
    1. refpanel_leg (list) - a list of snp (rsids, ref_alt) in the legend
    2. snp_beta (dict) - a dictionary that maps snps with beta
    3. snp_beta_info (list) - a list of (rsid, pos, n, ref_alt)
return:
    1. a dictionary mapping snp id with beta with sign fixed
    2. a list of filtered (snpid, pos, n, ref_alt)
"""
def filter_snps(refpanel_leg, snp_beta, snp_beta_info):
    # create a diction of snpid -> refalt
    snpid_refpanel_refalt = dict()
    for snpid,refalt in refpanel_leg:
        snpid_refpanel_refalt[snpid] = refalt
    
    # find snps to be removed
    filter_set = set()
    for i in xrange(len(snp_beta_info)):
        snpid = snp_beta_info[i][0]
        refalt = snp_beta_info[i][3]
        # check if snp is in reference panel
        if(snpid in snpid_refpanel_refalt):
            refpanel_refalt = snpid_refpanel_refalt[snpid]
            # snps with matching alleles
            if((refpanel_refalt == "AC" and (refalt == "TG" or refalt == "AC" or refalt == "TC" or refalt == "AG")) or
               (refpanel_refalt == "AG" and (refalt == "TC" or refalt == "AG" or refalt == "TG" or refalt == "AC")) or
               (refpanel_refalt == "CA" and (refalt == "GT" or refalt == "CA" or refalt == "GA" or refalt == "CT")) or
               (refpanel_refalt == "CT" and (refalt == "GA" or refalt == "CT" or refalt == "GT" or refalt == "CA")) or
               (refpanel_refalt == "TC" and (refalt == "AG" or refalt == "TC" or refalt == "AC" or refalt == "TG")) or
               (refpanel_refalt == "TG" and (refalt == "AC" or refalt == "TG" or refalt == "AG" or refalt == "TC")) or
               (refpanel_refalt == "GA" and (refalt == "CT" or refalt == "GA" or refalt == "CA" or refalt == "GT")) or
               (refpanel_refalt == "GT" and (refalt == "CA" or refalt == "GT" or refalt == "CT" or refalt == "GA"))):
                pass
            # snps with inverse alleles
            elif((refpanel_refalt == "AC" and (refalt == "GT" or refalt == "CA" or refalt == "CT" or refalt == "GA")) or
                 (refpanel_refalt == "AG" and (refalt == "CT" or refalt == "GA" or refalt == "GT" or refalt == "CA")) or
                 (refpanel_refalt == "CA" and (refalt == "TG" or refalt == "AC" or refalt == "AG" or refalt == "TC")) or
                 (refpanel_refalt == "CT" and (refalt == "AG" or refalt == "TC" or refalt == "TG" or refalt == "AC")) or
                 (refpanel_refalt == "TC" and (refalt == "GA" or refalt == "CT" or refalt == "CA" or refalt == "GT")) or
                 (refpanel_refalt == "TG" and (refalt == "CA" or refalt == "GT" or refalt == "GA" or refalt == "CT")) or
                 (refpanel_refalt == "GA" and (refalt == "TC" or refalt == "AG" or refalt == "AC" or refalt == "TG")) or
                 (refpanel_refalt == "GT" and (refalt == "AC" or refalt == "TG" or refalt == "TC" or refalt == "AG"))):
                snp_beta[snpid] = -1.0*snp_beta[snpid]
            # strand ambiguous snps
            elif(refalt == "AT" or refalt == "CG" or
                 refalt == "TA" or refalt == "GC"):
                filter_set.add(i)
            # snps with alleles not matching reference panel
            else:
                filter_set.add(i)
        # snp not found in reference panel
        else:
            filter_set.add(i)
    
    # create new snp_beta snp_beta_info
    snp_beta_filt = dict()
    snp_beta_info_filt = []
    for i in xrange(len(snp_beta_info)):
        if(i not in filter_set):
            snpid = snp_beta_info[i][0]
            snp_beta_filt[snpid] = snp_beta[snpid]
            snp_beta_info_filt.append(snp_beta_info[i])

    return snp_beta_filt,snp_beta_info_filt
