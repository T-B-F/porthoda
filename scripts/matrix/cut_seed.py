#!/usr/bin/env python
""" cut hmm file in multiple files
"""
import os, sys, argparse
import codecs

def get_cmd():
    """ read command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", action="store", dest="inputfile", help="file with multiple hmm models")
    parser.add_argument("-o", action="store", dest="dirout", help="directory to store separated hmm models")
    params = parser.parse_args()
    return params

def main():
    params = get_cmd()
    
    # create output dir if not existing
    if not os.path.isdir(params.dirout):
        os.makedirs(params.dirout)
        
    name = None
    lines = []
    data = dict()
    # read hmm files
    # the file from hhsuite database is often weirdly encoded latin-1 seems to solve it
    with open(params.inputfile, "r", encoding='latin-1') as input_handle:
        for line in input_handle:
            #print(line)
            if name:
                # the name if known, the line can be directly written in the correct output
                outf.write(line)
            else:
                # storing lines before knowing name of the models, ie name of output
                lines.append(line)
            if line.startswith("NAME "):
                # get name of the model and create the corresponding output file
                tmp = line.split()
                name = tmp[1].split(".")[0]
                idname = tmp[3]
                pathout = os.path.join(params.dirout, name+".hmm")
                outf = open(pathout, "w")
                for line in lines:
                    # write down previously stored lines
                    outf.write(line)
            elif line.startswith("//"):
                name = None
                lines = []
                outf.close()
        
    sys.exit(0)

if __name__ == "__main__":
    main()