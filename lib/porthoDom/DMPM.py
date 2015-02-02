#!/usr/bin/env python
""" DMPM module 

implement a DMPM class to make the interface with the similarity matrix
"""

import sys
import struct


__author__ = "Carste Kemena, Tristan Bitard Feildel"
__email__ = "t.bitard.feildel@uni-muenster.de"
__institute__ = "Insitute for Evolution and Biodiversity"
__lab__ = "Evolutionary Bioinformatics"
__vesion__ = "1.0"


class DMPM:
    """
    Domain Match Probability Matrix
    Read a domain similarity matrix stored into a binary parse column format 
    row and column indexes of the matrix are the id of the domains models
    Each matrix are created according to a cutoff, value below this cutoff are 
    not stored inside the matrix and a -1 similarity score will be return if 
    queried
    """
    def __init__(self, path):
        """init function of DMP matrix 
        
        Parameter
        ---------
            path : string
                input file, path to the binary matrix
        """
        self.n_domains = 0
        self.n_vals = 0
        with open(path, "rb") as in_f:
            self.name = struct.unpack('10s', in_f.read(10))[0]
            self.n_domains = struct.unpack('i', in_f.read(4))[0]
            self.n_vals = struct.unpack('i', in_f.read(4))[0]
            ids = struct.unpack('i'*self.n_domains, in_f.read(4*self.n_domains))
            self.ids = {}
            for i in range(0, len(ids)):
                self.ids[ids[i]] = i
            self.row_ids = struct.unpack('i'*self.n_domains, 
                                            in_f.read(4*self.n_domains))
            self.col_ids = struct.unpack('i'*self.n_vals, 
                                            in_f.read(4*self.n_vals))
            self.values = list(struct.unpack('h'*self.n_vals, 
                                            in_f.read(2*self.n_vals)))
    
    def get_val(self, idi, idj):
        """Return a value of the matrix
        
        Parameters
        ----------
            idi : int 
                domain identification number
            idj : int 
                domain identification number
        Return
        ------
            value : int
                the similarity value, or -1 if not found 
                (ie similarity < cutoff matrix)
        
        """
        id1 = -1
        id2 = -1
        if idi == idj:
            raise ValueError("both id are indentical")
        if idi < idj:
            id1 = self.ids[idi]
            id2 = self.ids[idj]
        else:
            id1 = self.ids[idj]
            id2 = self.ids[idi]
        id_col = self.row_ids[id1]
        if id_col != (self.n_domains-1):
            end = self.row_ids[id1+1]
        else:
            end = self.n_vals
        for i in range(id_col, end):
            if self.col_ids[i] == id2:
                return self.values[i]
        return -1

    def set_val(self, idi, idj, val):
        """Change a value inside the matrix
        
        Parameters
        ----------
            idi : int 
                domain identification number
            idj : int 
                domain identification number
            val : int 
                new value
        Return
        ------
            value : int
                0 if the value has been changed, -1 otherwise
        
        """
        id1 = -1
        id2 = -1
        if idi < idj:
            id1 = self.ids[idi]
            id2 = self.ids[idj]
        else:
            id1 = self.ids[idj]
            id2 = self.ids[idi]

        id_col = self.row_ids[id1]
        if id_col != (self.n_domains-1):
            end = self.row_ids[id1+1]
        else:
            end = self.n_vals
        for i in range(id_col, end):
            if self.col_ids[i] == id2:
                self.values[i] = val
                return 0
        return -1

def main():
    """ Test function
    """
    pathmat = sys.argv[1]
    val1 = int(sys.argv[2])
    val2 = int(sys.argv[3])
    val = DMPM(pathmat)
    print val.get_val(val1, val2)


if __name__ == "__main__":
    main()

