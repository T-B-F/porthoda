# -*- coding: utf-8 -*-
#!/usr/bin/env python
# I/O functions for proteinorthoDA.py
# Copyright (C) 2013  Institute for  Evolution and Biodiversity
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

__author__ = "Tristan Bitard Feildel"
__email__ = "t.bitard.feildel@uni-muenster.de"
__institute__ = "Insitute for Evolution and Biodiversity"
__lab__ = "Evolutionary Bioinformatics"
__vesion__ = "1.0"

import sys, pickle
import networkx as nx 

__all__ = [ "extractPFidDA", "extract_PFid_to_xdom", "read_multifasta",
"read_porthoparams", "write_results_daclusters", "write_results_porthodaclusters",
"read_dadone" ]


def read_dadone( path, luniq_da ) :
    """ Read the already computed DA similarity and create the graph, dict and 
    approriated list
    
    Parameters
    ----------
    path : string
        da similarity file
    luniq_da : list
        the list of unique domain arrangement
        
    Returns
    -------
    GDA : nx.Graph( )
        DA similarity graph
    da_similarity : dict 
        pairwise similarity in graph format
     missing : list
        domain arrangement in luniq_da but not in the graph
    """
    h = open( path ) 
    da_similarity = pickle.load( h ) 
    h.close( )
    GDA = nx.Graph( )
    for k1 in da_similarity :
        for k2 in da_similarity[k1] :
            GDA.add_edge( k1, k2, weight=da_similarity[k1][k2] )
    missing = [ ]
    for da in luniq_da :
        if not da in GDA :
            missing.append( da )
    return GDA, da_similarity, missing
    
    
    

def suppress_repeat( dom_list ) :
    """Suppress direct repeat in a list
    suppress_repeat( [ a,a,b,b,c] ) -> [ a,b,c ]
    
    Parameter
    ---------
        dom_list : the list with repeat
        
    Return
    ------
        new_list : the list without repeat
    
    """
    new_list = [dom_list[0]]    
    for i in range(1, len(dom_list)):
        if dom_list[i] != new_list[-1] :
            new_list.append( dom_list[i] )
    return new_list

def extractPFidDA( pfamscanfile , replace_by_clan = False, mask_repeat=False) :
    """Extract Pfam id for domain arrangement
    
    Parameters
    ----------
        pfamscanfile : the path to the pfam output file
        replace_by_clan : if yes use Clan name if avaliable instead of Pfam id
        
    Return
    ------
        dprot a dictionary [ protein ] [ domain_arrangement_1, ... ]
    """
    h = open( pfamscanfile ) 
    line = h.readline( )
    dprot = { }
    while line != "" :
        if line[0] == "\n" or line[0] == "#" : 
            line = h.readline( )
            continue 
        tmp = line.split( )
        dom = tmp[5].split(".")[0]
        if replace_by_clan :
            if tmp[14] != "No_clan" :
                dom = tmp[14]
        ret = dprot.setdefault( tmp[0], [] ).append( dom )
        line = h.readline( )
    if mask_repeat :
        for prot in dprot :
            doms = dprot[prot]
            new_doms = suppress_repeat( doms )
            dprot[ prot ] = new_doms 
    return dprot


def extract_PFid_to_xdom( pfamscanfile , replace_by_clan = False) :
    """Convert a pfam output file to an xdom file
    
    Parameters
    ----------
        pfamscanfile : the path to the pfam output file
        replace_by_clan : if yes use Clan name if avaliable instead of Pfam id
    
    Returns
    -------
        dprot : a dictionary [ protein ] [ domain_arrangement_1, ... ]
    """
    # Nested domains are removed, this file will be used by rads/rampage and 
    # currently the program doesn't handle nested domains
    h = open( pfamscanfile ) 
    line = h.readline( )
    dprot = { }
    while line != "" :
        if line[0] == "\n" or line[0] == "#" : 
            line = h.readline( )
            continue 
        tmp = line.split( )
        dom = tmp[5].split(".")[0]
        start = int(tmp[1])
        end = int( tmp[2] )
        if replace_by_clan :
            if tmp[14] != "No_clan" :
                dom = tmp[14]
        ret = dprot.setdefault( tmp[0], [] ).append( (dom,start,stop) )
        line = h.readline( )
    return dprot

def read_multifasta( filename ) :
    """Read the a multifasta file
    
    Parameters
    ----------
        filename : path to the multifasta file
    Return 
    ------
        d : a dictionary associating d[ fasta header ] = fasta sequence
    """
    d = { }
    head = None 
    f = open(filename)
    for line in f :
        if line[0] == ">" :
            if head != None :
                d[head] = seq
            head = line[1:].split()[0] # fasta header name is the first entry before a space
            seq = ""
        else :
            seq += line.strip()
    if line[0] != ">" and seq != "":
        d[head] = seq
    else :
        raise ValueError,"Expected an amino acid sequence in last line"
    f.close( )
    return d 

def read_porthoparams( pathfile, porthopath ) :
    """ Read a file with proteinortho.pl user specific parameter
    If no file the default parameters are used otherwise
    
    Parameters
    ----------
        pathfile : path to the file containing the parameters
        porthopath : path to proteinortho binary

    Return 
    ------
        portho_params : list of proteinortho params
    """
    portho_params = [ ]    
    list_params = ["-e","-p","-id","-cov","-conn","-m","-pairs","-singles",
    "-selfblast","-unambiguous","-a","-f","-ff"]
    need_value = [1,1,1,1,1,1,0,0,0,0,1,0,0]
    incompatible = ["-dir","-remove","-batch","-blastonly","-blastdone","-log",
    "-o","-verbose","-debug" ]
    h = open(pathfile)
    f = h.readlines()
    h.close( )
    for line in f :
        param = line[:-1]
        tmp = param.split("=")
        if tmp[0] in incompatible :
            print >>sys.stderr, "Error : Incompatible proteinortho parameters (%s) in %s"%(param,pathfile)
            sys.exit(1)
        if tmp[0] not in list_params :
            print >>sys.stderr, "Error : Unknown parameters %s for "+porthopath+" in %s"%(param,pathfile)
            sys.exit( 1 )
        if "=" in param and tmp[1] == "" :
            print >>sys.stderr, "Error : Parameters (%s) in %s need a parameter value"%(param,pathfile)
            print >>sys.stderr, "        Check '"+porthopath+" -h' for more information"
            sys.exit(1)
        i = list_params.index( tmp[0] )
        if need_value[i] == 1 and len( tmp ) == 1 :
            print >>sys.stderr, "Error : Parameters (%s) in %s need a parameter value"%(param,pathfile)
            print >>sys.stderr, "        Check '"+porthopath+" -h' for more information"
            sys.exit(1)
        portho_params.append( params )
    return portho_params
                                                                                                                                                                                                                                                                       
    
def write_results_daclusters( clusters, uniq_da, name_proteomes, param ) :
    """Output results of DA clustering in a proteinortho fileformat way
    
    Parameters
    ----------
        clusters : list[list]
            a list of list containing da clustered
        uniq_da : dict
            dictionary mapping domain arrangement and protein/specie
        name_proteome : dict
            mapping between internal names of proteomes and external names
        param : argument parser object
            store options relative to the program
    """
    lsp = name_proteomes.keys( )
    dresults = { }
    portho_info = { }
    tot_prot = 0
        
    for cnt, compo in enumerate( clusters ):       
        dresults[cnt] = { }
        for da in compo :
            spprots = uniq_da[ da ] # uniq_da is the dictionary associating a domain arrangement to a list of species/proteins tuple
            for spprot in spprots :  
                sp,prot = spprot.split(";")
                #print spprot
                dresults[cnt].setdefault(sp,[]).append( prot )                    
                ret = portho_info.setdefault( int(sp), 0 )
                portho_info[ int(sp) ] += 1 
                tot_prot += 1 
            
    # create one single file
    sortie = file( param.output, "w" )
    sortie.write("#species\tproteins\talg.-conn.")
    for sp in lsp :
        sortie.write("\t"+name_proteomes[sp]) 
        
    sortie.write( "\n" )
    sortie.write("#"+str(len(portho_info.keys()))+"\t"+str(tot_prot)+"\t-")
    
    for sp in lsp :
        if portho_info.has_key( sp ) : # no key mean no protein of this specie found in a cluster
            sortie.write("\t"+str(portho_info[sp]))
        else :
            sortie.write("\t0")
    sortie.write("\n")
        
    for cnt in dresults :
        nbspecies = len(dresults[cnt].keys())
        nbprots = sum( [ len(val) for val in dresults[cnt].values() ] )
        new_line = str(nbspecies) + "\t"+ str(nbprots) + "\tNan"
        for sp in lsp :
            if dresults[cnt].has_key( str(sp) ) :
                prots = ",".join( dresults[cnt][str(sp)] )
                new_line += "\t" + prots
            else :
                new_line += "\t*"
        sortie.write( new_line + "\n" )
    sortie.write("#cutoff:{}\torder:{}\trepeat:{}\tepsilon:{}\tminpts:{}\n"\
                 .format(param.cutoff, param.order,param.maskrepeat,
                         param.epsilon,param.minpts) )
    sortie.close( )
    

def write_results_porthodaclusters( name_proteomes, portho_info, dresults, footer, param  ) :
    """Output results of DA + proteinortho clustering in a proteinortho fileformat way
    
    Parameters
    ----------
        name_proteome : dict
            mapping between internal names of proteomes and external names
        portho_info : dict
            store information relative to proteinortho (number of specie, number of prot/per specie)
        dresults : dict
            result storing result of proteinortho runs
        footer : string
            option and version of proteinortho program
        param : argument parser object
            store options relative to the program
    """
    tot_prot = 0
    lsp = sorted( name_proteomes.keys() )
    for sp in portho_info.keys( ) :
        tot_prot += len( portho_info[sp].keys( ) )

    # create one single file
    sortie = file( param.output, "w" )
    sortie.write("#species\tproteins\talg.-conn.")
    
    for sp in lsp :
        sortie.write("\t"+name_proteomes[sp])
        
    sortie.write("\n")
    sortie.write("#"+str(len(lsp))+"\t"+str(tot_prot)+"\t-")
    
    for sp in lsp :
        if portho_info.has_key( sp ) : # no key mean no protein of this specie found in a cluster
            sortie.write("\t"+str(len(portho_info[sp].keys())))
        else :
            sortie.write("\t0")
            
    sortie.write("\n")
    
    for fam in dresults :
        nbspecies = dresults[ fam ][ 0 ]
        nbprot = dresults[ fam ][ 1 ]
        nbconn = dresults[ fam ][ 2 ]
        prots = dresults[ fam ][ 3 ]
        new_line = nbspecies + "\t"+nbprot + "\t" + nbconn
        species,lprots = zip(*prots)
        
        for sp in lsp :
            if sp in species :
                ind = species.index( sp )
                new_line += "\t" + lprots[ind]
                
            else :
                new_line += "\t*"
                
        sortie.write( new_line + "\n" )
        
    sortie.write(footer.strip()+" cutoff:{}\torder:{}\trepeat:{}\tepsilon:{}\tminpts:{}\n"\
                 .format(param.cutoff, param.order,param.maskrepeat,
                         param.epsilon,param.minpts) )
    sortie.close( )