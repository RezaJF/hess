#!/usr/bin/python
# (c) 2017-2022 Huwenbo Shi

import argparse, math, sys, logging, time

# main function
def main():
    
    # get command line input
    args = get_command_line()

    # if shared sample size is 0 output 0
    if args.ns == 0.0:
        print 0.0

    # parse cross-trait LDSC log
    phenocor = None
    log_f = open(getattr(args, 'ldsc-log'), 'r')
    counter = 0
    for line in log_f:
        line = line.strip()
        if counter > 0: counter += 1
        if line.strip() == 'Genetic Covariance':
            counter += 1
        if counter == 5:
            intercept = float(line.split()[1])
            phenocor = intercept * math.sqrt(args.n1*args.n2) / args.ns
            break
    log_f.close()

    # print the results
    if phenocor != None:
        print phenocor

# get command line
def get_command_line():
    
    parser = argparse.ArgumentParser(description='Estimate phenotypic '\
        'correlation from cross-trait LDSC intercept')

    parser.add_argument('--ldsc-log', dest='ldsc-log', type=str,
        help='Cross-trait LDSC log file', required=True)
    
    parser.add_argument('--n1', dest='n1', type=float,
        help='Sample size for GWAS 1', required=True)

    parser.add_argument('--n2', dest='n2', type=float,
        help='Sample size for GWAS 2', required=True)

    parser.add_argument('--ns', dest='ns', type=float,
        help='Number of shared samples', required=True)

    args = parser.parse_args()
    
    return args

if(__name__ == '__main__'):
    main()
