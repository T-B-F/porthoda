#!/usr/bin/env python
""" read each file and create the appropriated matrix in CRS format
CRS format code adapted from Carsten Kemena C++ code
"""
from math import log, floor
import argparse
import os, sys
import numpy as np
import struct

def getData(line, name):
    """ get values from hhsearch similarity comparison
    """
    tmp = line.split()
    # depending on input file, results of hhsearch have some time an intermetiate 
    # column after the domain name
    if line[10:30].split() == []:
        offset = 4
    else :
        offset = 5
    # get the value of the line
    if name == "LOG-PVA":
        val = float(tmp[offset])
    elif name == "S-AASS":
        val = float(tmp[offset+1])
    elif name == "PROBAB":
        val = float(tmp[offset+2])
    elif name == "SCORE":
        val = float(tmp[offset+3])
    elif name == "LOG-EVAL":
        val = float(tmp[offset+4])
    else:
        print("Error, unknow name to extract {}".format(name), file=sys.stderr)
        sys.exit(1)
    return int(floor(val+ 0.5))

def prefix_file(path):
    """ get name of a file without extension out of a path
    """
    return os.path.splitext(os.path.basename(path))[0]

def get_cmd():
    """ read command line arguments
    """
    parser = argparse.ArgumentParser( )
    parser.add_argument("-i", action="store", dest="listres", nargs="+", 
                        help="list of hhsearch score files")
    parser.add_argument("-n", action="store", dest="name", 
                        help="matrix name")
    parser.add_argument("-s", action="store", dest="scorename", 
                        help="score name to keep", choices=["LOG-PVA", 
                            "S-AASS","PROBAB", "SCORE", "LOG-EVAL"])
    parser.add_argument("-t", action="store", dest="threshold", type=float,
                        help="threshold applied to the score", default=1)
    parser.add_argument("-c", action="store", dest="cutoff", 
                        help="cutoff evalue", type=float, default=0.1)
    parser.add_argument("-o", action="store", dest="outmatrix", 
                        help="output matrix")
    params = parser.parse_args()
    assert 0 < len(params.name)  <= 10, "Error length of matrix name should be comprised between 0 and 11"
    return params

def main():
    # read arguments
    params = get_cmd()
    
    # read mapping
    mapping_name2id = read_mapping(params.mapping)

    # pfam num id
    # TODO need to be change for other databases
    names, inputfiles = list(zip(*sorted([(os.path.basename(f).split(".")[0], f) for f in params.listres])))
    # mapping pfamid, readl idx (some value can be missing in domain numbering
    dnames = dict(zip(names, range(len(names))))
       
    vals, allids, col_ids, row_ids = [], [], [], []
    for row_num, inputf in enumerate(inputfiles):
        name = names[row_num]
        # read data
        allids.append(int(name[2:]))
        row_ids.append(len(col_ids))
        with open(inputf) as inf:
            found = False
            for line in inf :
                tmp = line.split()
                if found:
                    target = tmp[0]
                    col_num = dnames[target]
                    if col_num >= row_num:
                        val = getData(line, params.scorename)
                        evalue = getData(line, "LOG-EVAL")/(-1.443) # evalue formula in HHsearch output -1.443 * hit.logEval                     
                        if evalue < log(params.cutoff) and val > params.threshold:
                            vals.append(val)
                            col_ids.append(col_num)
                elif line.startswith("TARGET"):
                    found = True
   
    ndomains = len(names)
    nvals = len(col_ids)
    
    with open(params.outmatrix, "wb") as outf:
        outf.write(struct.pack("<10s", str.encode(params.name, "utf-8")))
        outf.write(struct.pack("<i", ndomains))
        outf.write(struct.pack("<i", nvals))
        outf.write(struct.pack("<"+"i"*ndomains, *allids)) 
        outf.write(struct.pack("<"+"i"*ndomains, *row_ids))
        outf.write(struct.pack("<"+"i"*nvals,    *col_ids))
        outf.write(struct.pack("<"+"h"*nvals,    *vals))
        
    sys.exit( 0 )
    
if __name__ == "__main__":
    main()