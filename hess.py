#!/usr/bin/python
# (c) 2017-2022 Huwenbo Shi

import numpy as np
import pandas as pd
import argparse, math, sys, logging, time
from pysnptools.snpreader import Bed
from src.estimation import *

title_str = """\
@----------------------------------------------------------@
       |         HESS       |      v0.5      |    9/October/2017  |
       |----------------------------------------------------------|
       |  (C) 2017 Huwenbo Shi, GNU General Public License, v3    |
       |----------------------------------------------------------|
       |  For documentation, citation & bug-report instructions:  |
       |   http://bogdan.bioinformatics.ucla.edu/software/hess/   |
       @----------------------------------------------------------@\
"""

# main function
def main():

    # get command line argument and initialize log
    args = get_command_line()
    init_log(args)

    # get analysis code
    analysis, argmap = check_command_line(args)

    # run the corresponding analysis
    if analysis == 'local_hsqg_step1':
        local_hsqg_step1(argmap['bfile'], argmap['local-hsqg'],
            argmap['partition'], argmap['chrom'], argmap['min-maf'],
            argmap['out'])
    
    elif analysis == 'local_hsqg_step2':
        local_hsqg_step2(argmap['prefix'], argmap['max-num-eig'],
            argmap['min-eigval'], argmap['reinflate-lambda-gc'],
            argmap['gwse-thres'], argmap['tot-hsqg'], argmap['out'])
    
    elif analysis == 'local_rhog_step1':
        local_rhog_step1(argmap['bfile'], argmap['local-rhog'],
            argmap['partition'], argmap['chrom'], argmap['min-maf'],
            argmap['out'])
    
    elif analysis == 'local_rhog_step2':
        local_rhog_step2(argmap['prefix'], argmap['local-hsqg-est'],
            argmap['max-num-eig'], argmap['min-eigval'],
            argmap['reinflate-lambda-gc'], argmap['gwse-thres'],
            argmap['pheno-cor'], argmap['num-shared'], argmap['out'])

    else:
        logging.error('Invalid command line option.')

    # end the log
    end_log()

# check command line
def check_command_line(args):

    # parse command line
    argmap = dict()
    for arg in vars(args):
        argmap[arg] = getattr(args, arg)

    # step 1 for local SNP-heritability
    if(argmap['chrom'] != None and argmap['partition'] != None and
       argmap['local-hsqg'] != None and argmap['bfile'] != None):
        return ('local_hsqg_step1', argmap)

    # step 2 for local SNP-heritability
    if(argmap['prefix'] != None and argmap['out'] != None and
       argmap['pheno-cor'] == None and argmap['num-shared'] == None and
       argmap['local-hsqg-est'] == None):
        return ('local_hsqg_step2', argmap)

    # step 1 for local genetic covariance
    if(argmap['chrom'] != None and argmap['partition'] != None and
       argmap['local-rhog'] != None and argmap['bfile'] != None):
        return ('local_rhog_step1', argmap)

    # step 2 for local genetic covariance
    if(argmap['prefix'] != None and argmap['out'] != None and
       argmap['pheno-cor'] != None and argmap['num-shared'] != None and
       argmap['local-hsqg-est'] != None):
        return ('local_rhog_step2', argmap)

    # command line option not recognized
    return ('invalid', argmap)

# initialize log
def init_log(args):

    # get log file name
    log_file_name = args.out
    if args.chrom != None:
        log_file_name = log_file_name + '_chr' + args.chrom
    log_file_name += '.log'

    # create the log file
    log_format = '[%(levelname)s] %(message)s'
    logging.basicConfig(filename=log_file_name, filemode="w",
        level=logging.DEBUG, format=log_format)

    # add stderr as a stream handler
    stderr_handler = logging.StreamHandler(sys.stderr)
    formatter = logging.Formatter(log_format)
    stderr_handler.setFormatter(formatter)
    logging.getLogger().addHandler(stderr_handler)

    # log time and command issued
    logging.info(title_str)
    specified = set([val for val in sys.argv if val[0] == '-'])
    cmd_str = sys.argv[0] + ' \\\n'
    for arg in vars(args):
        if ('--' + arg) in specified:
            param = getattr(args, arg)
            if type(param) == list:
                param = ' '.join([str(p) for p in param])
            elif type(param) == bool:
                param = ''
            cmd_str += '        --{} {} \\\n'.format(arg, param)
    cmd_str = cmd_str.strip()[0:len(cmd_str)-3]
    cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    logging.info('Command started at: %s' % cur_time)
    logging.info('Command issued:\n    {}'.format(cmd_str))

# end the log
def end_log():
    cur_time = time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
    logging.info('Command finished at: %s' % cur_time)

# get command line input
def get_command_line():
    
    # create the help document
    parser = argparse.ArgumentParser(description='Estimate local '
        'SNP-heritability and genetic covariance from GWAS summary data')

    ########## step 1 ##########

    parser.add_argument('--bfile', dest='bfile', type=str, required=False,
        default=None, help='Reference panel file in PLINK file format')

    parser.add_argument('--chrom', dest='chrom', type=str, required=False,
        help='Specifies the chromosome number')
    
    parser.add_argument('--partition', dest='partition', type=str,
        help='A bed format file specifying how the genome is partitioned',
        required=False)
    
    parser.add_argument('--local-hsqg', dest='local-hsqg', type=str,
        help='Z-score file for local SNP-heritability estimation',
        required=False)

    parser.add_argument('--local-rhog', dest='local-rhog', type=str, nargs=2,
        help='Z-score files for local genetic covariance estimation',
        required=False)

    ########## step 2 ##########

    # local SNP-heritability

    parser.add_argument('--tot-hsqg', dest='tot-hsqg', nargs=2, type=float,
        help='Total trait SNP heritability and standard error', required=False)
    
    # local genetic covariance

    parser.add_argument('--num-shared', dest='num-shared', type=float,
        help='Number of shared samples', required=False)
    
    parser.add_argument('--pheno-cor', dest='pheno-cor', type=float,
        help='Phenotype correlation', required=False)
    
    parser.add_argument('--local-hsqg-est', dest='local-hsqg-est', type=str,
        nargs=2, default=None, help='Local SNP-heritability estimates')
    
    ########## other arguments ##########
    
    # prefixe used for step 1
    parser.add_argument('--prefix', dest='prefix', type=str,
        help='Prefix used for step 1', required=False)

    # maximum number of eigenvectors use in the truncated-SVD regularization
    # of the LD matrix used
    parser.add_argument('--max-num-eig', dest='max-num-eig', type=int, 
        default=50, help='Maximum number of eigenvectors to use in the'
        'truncated-SVD regularization of the LD matrix (default 50)')

    # eigenvalue threshold
    parser.add_argument('--min-eigval', dest='min-eigval', type=float,
        default=1.0, help='Eigenvalue threshold (default 1.0)')

    # maf threshold
    parser.add_argument('--min-maf', dest='min-maf', type=float,
        default=0.05, help='Minimum minor allele frequency')

    # refinlate the estimates by the lambda gc
    parser.add_argument('--reinflate-lambda-gc', nargs='*',
        dest='reinflate-lambda-gc', type=float, required=False, default=None,
        help='Reinflate the summary stats with by the genomic control factor')

    # a threshold on the sensitivity of the genome-wide estimates
    parser.add_argument('--gwse-thres', dest='gwse-thres', type=float,
        default=None, help='Sensitivity threshold on total SNP-heritability'
        ' estimates (default is 2.0 if --tot-hsqg is not specified and 0.5'
        ' otherwise')

    # specify output file
    parser.add_argument('--out', dest='out', type=str,
        help='Output file name', required=True)
    
    return parser.parse_args()

if(__name__ == '__main__'):
    main()
