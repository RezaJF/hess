#!/usr/bin/python

from optparse import OptionParser
import numpy as np, numpy.linalg
import os, sys, math

eps = 10**-8

# main function
def main():
   
    # get command line
    params = get_command_line()
    prefix = params["prefix"]
    thres = params["thres"]
    out_file = params["out_file"]
    gc = params["gc"]

    # load step1
    locus_info,all_eig,all_prj = load_step1(prefix)
    
    # get h2g
    all_h2g,raw_est = get_local_h2g(locus_info, all_eig, all_prj, thres, gc)

    # get variance estimates
    all_var = get_var_est(locus_info, all_h2g)

    # write output
    out_file = open(out_file, 'w')
    for i in xrange(len(locus_info)):
        line = '%s\t%s\t%s\t%s\t%s\t%d\t%s\t%.10f\t%.10f\t%.12f' % (
                locus_info[i][0], locus_info[i][1], locus_info[i][2],
                locus_info[i][3], locus_info[i][4], raw_est[i][1],
                locus_info[i][5], raw_est[i][0], all_h2g[i,0], all_var[i])
        out_file.write(line+'\n')
    out_file.close()

# compute heritability
def get_local_h2g(locus_info, all_eig, all_prj, thres, gc):
    raw_est = []
    for i in xrange(len(locus_info)):
        n = float(locus_info[i][5])
        k = min(thres, all_eig[i].size)
        tmp = np.divide(all_prj[i][0,0:k], all_eig[i][0,0:k]+eps)
        raw_est.append((np.sum(tmp), float(k)))
    num_win = len(locus_info)
    est = np.matrix(np.zeros((num_win, 1)))
    for i in xrange(num_win):
        n = float(locus_info[i][5])
        est[i,0] = n*raw_est[i][0]*gc-raw_est[i][1]
        est[i,0] /= n-raw_est[i][1]+eps
    return est,raw_est

# est var
def get_var_est(locus_info, all_h2g):
    num_win = len(locus_info)
    tot = np.sum(all_h2g)
    var_est = []
    for i in xrange(num_win):
        h2g = all_h2g[i,0]
        n = float(locus_info[i][5])
        p = float(locus_info[i][4])
        var = 1.0
        var = var*((n/(n-p+eps))**2)
        var = var*(2*p*((1-h2g)/(n+eps))+4*h2g)
        var = var*((1-h2g)/(n+eps))
        var_est.append(var)
    return var_est

# load step 1
def load_step1(prefix):
    
    # load info
    locus_info = []
    for i in xrange(1,23):
        fnm = '%s_chr%d.info' % (prefix, i)
        if(not os.path.exists(fnm)):
            continue
        fnm = open(fnm)
        for line in fnm:
            line = line.strip()
            cols = line.split()
            locus_info.append((i,cols[0],cols[1],cols[2],cols[3],cols[4]))
        fnm.close()

    # eigs
    all_eig = []
    for i in xrange(1,23):
        fnm = '%s_chr%d.eig' % (prefix, i)
        if(not os.path.exists(fnm)):
            continue
        fnm = open(fnm)
        for line in fnm:
            line = line.strip()
            cols = line.split()
            tmp = np.matrix([float(cols[i]) for i in range(len(cols))])
            all_eig.append(tmp)
        fnm.close()
    
    # prjsq
    all_prj = []
    for i in xrange(1,23):
        fnm = '%s_chr%d.prjsq' % (prefix, i)
        if(not os.path.exists(fnm)):
            continue
        fnm = open(fnm)
        for line in fnm:
            line = line.strip()
            cols = line.split()
            tmp = np.matrix([float(cols[i]) for i in range(len(cols))])
            all_prj.append(tmp)
        fnm.close()
    
    return locus_info,all_eig,all_prj

# get command line
def get_command_line():
    
    # get command line
    parser = OptionParser(add_help_option=False)
    parser.add_option("-p", "--prefix", dest="prefix")
    parser.add_option("-t", "--threshold", type="int",
        dest="thres", default=30)
    parser.add_option("-o", "--out_file", dest="out_file")
    parser.add_option("-g", "--gc", dest="gc", type="float",
        default=1.0)
    (options, args) = parser.parse_args()
    prefix = options.prefix
    thres = options.thres
    out_file = options.out_file
    gc = options.gc

    # check command line
    if(prefix == None or thres == None or out_file == None):
        sys.stderr.write("Usage:\n")
        sys.stderr.write("\tUse -p to specify prefix\n")
        sys.stderr.write("\tUse -t to specify threshold\n")
        sys.stderr.write("\tUse -o to specify output file\n")
        sys.exit()
    
    # create param dictionary
    params = dict()
    params["prefix"] = prefix
    params["thres"] = int(thres)
    params["out_file"] = out_file
    params["gc"] = float(gc)

    return params

if(__name__ == '__main__'):
    main()
