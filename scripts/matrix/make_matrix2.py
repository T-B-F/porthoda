#!/usr/bin/env python
""" read each file and create a full matrix
use mat2bincrs to transform it
"""
from math import log
import os, sys
import numpy as np
import argparse

def prefix_file(path):
    """ get name of a file without extension out of a path
    """
    return os.path.splitext(os.path.basename(path))[0]

def getData(line, name):
    """ get value 
    """
    tmp = line.split()
    # depending on input file, results of hhsearch have some time an intermetiate 
    # column after the domain name
    if line[10:30].split() == [] :
        offset = 4
    else :
        offset = 5
    if name == "LOG-PVA" :
        val = float( tmp[offset] )
    elif name == "S-AASS" :
        val = float( tmp[offset+1]  )
    elif name == "PROBAB" :
        val = float( tmp[offset+2]  )
    elif name == "SCORE" :
        val = float( tmp[offset+3]  )
    elif name == "LOG-EVAL" :
        val = float( tmp[offset+4]  )
    else :
        print >>sys.stderr,"Unknow name to extract "
        sys.exit(1)
    return val

def read_mapping(path):
    """ read mapping between id and names
    """
    data = dict()
    with open(path) as inf:
        for line in inf:
            tmp = line.split()
            # key is the name, value is the id, ex( data["7tm_1"] = "PF00001"
            data[tmp[1]] = tmp[0]
    return data

def get_cmd():
    """ read command line parameters
    """
    parser = argparse.ArgumentParser( )
    parser.add_argument("-i", action="store", dest="listres", help="directory to result", nargs="+")
    parser.add_argument("-d", action="store", dest="dataname", help="data to get")
    parser.add_argument("-m", action="store", dest="mapping",  help="mapping between names and id")
    parser.add_argument("-c", action="store", dest="cutoff_evalue", help="cutoff evalue", type=float, default=10**(-3))
    parser.add_argument("-o", action="store", dest="outmatrix", help="output matrix")
    params = parser.parse_args( )
    return params

def main():
    # read command line
    params = get_cmd()
    
    # read mapping
    mapping_name2id = read_mapping(params.mapping)
    
    # get result names
    list_res = [(prefix_file(f), f) for f in params.listres]
    names, inputfiles = zip(*sorted(list_res))
    dnames = dict(zip(names, range(len(names)))
       
    mat = np.zeros((len(names), len(names)), dtype=np.float32)
    for i, inputf in enumerate(inputfiles):
        name = names[i]
        # read data 
        with open(inputf) as inf:
            found = False
            for line in f:
                if found:
                    tmp = line.split()
                    target = tmp[0][2:]
                    j = dnames[mapping_name2id[target]]
                    val = getData(line, params.dataname)
                    evalue = getData(line, "LOG-EVAL")/(-1.443) # evalue formula in HHsearch output -1.443 * hit.logEval 
                    if evalue < log(params.cutoff_evalue):
                        mat[i,j] = val
                elif line.startswith("TARGET"):
                    found = True

    with open(params.outmatrix, "w") as outf:
        line = ""
        for name in names:
            line += name+" "
        outf.write(line[:-1] + "\n")

        for i in range(mat.shape[0]):
            line = "" 
            for j in range(mat.shape[1]):
                line += str(mat[i,j])+" "
            outf.write(line[:-1] + "\n")

    sys.exit(0)
    
if __name__ == "__main__":
    main()
    