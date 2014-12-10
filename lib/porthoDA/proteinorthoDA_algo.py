# -*- coding: utf-8 -*-
#!/usr/bin/env python
# proteinorthoDA is a wrapper around proteinortho.pl script using a 
# pre-clustering step based on domain arrangement similarity to speed up the 
# blast all against all step used by proteinorhto.
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
""" proteinorthoDA_algo holds the main algorithms functions for clustering
and analysis
"""

import numpy as np
import networkx as nx

from proteinorthoDA_err import DASimilarityError, ProteinorthoError
from proteinorthoDA_util import timestamp, error_clean
from proteinorthoDA_worker import  subprocess_threaded_blastonly
from proteinorthoDA_worker import  subprocess_threaded_sim

import Pycluster
from Optics import Optics


__author__ = "Tristan Bitard Feildel"
__email__ = "t.bitard.feildel@uni-muenster.de"
__institute__ = "Insitute for Evolution and Biodiversity"
__lab__ = "Evolutionary Bioinformatics"
__vesion__ = "0.5"

__all__ = ["run_proteinortho_blast", "compute_similarity", "cluster_domains" ]

def run_proteinortho_blast ( clusters, uniq_da, ddfasta, params) :
    """Execute the blast part of proteinortho.pl 

    run the blast part of proteinortho either followed by a blastdone (ie 
    full analysis without any -blastonly or -blastdone parameter) or specificaly
    with the -blastonly argument (-blastdone should be call in a other run)
    
    Parameters
    ----------
    clusters : list 
        the list of clusters found in the graph
    uniq_da : dict 
        the dictionary linking da to a list of proteins
    ddfasta : dict 
        the dictionary containing the fasta sequence of all species and prot
        an empty object with specific parameters
    params : argparse.Namespace
        the dictionary containing the fasta sequence of all species and prot
        an empty object with specific parameters
    
    Returns
    -------
    listsub : list 
        the list of all submitted clusters ( for blastdone part )
    """
    listsub = [ ]
    dpid = { }
    lp = [ ]
    pool_blastonly = multiprocessing.Pool( processes = params.nb_job )
    dicres = { }
    
    for cnt, compo in enumerate( clusters ):
        subfamily = [ ]
        
        for da in compo :
            spprots = uniq_da[ da ] # uniq_da is the dictionary associating a domain arrangement to a list of species/proteins tuple
            
            for spprot in spprots :  
                subfamily.append( spprot )
                
        dportho = os.path.join( params.portho_dir , "sub" + str( cnt ) )
        dres = os.path.join( params.res_dir , "sub" + str( cnt ) )
        
        # create a sub directory per components
        if os.path.isdir ( dportho ) == False :
            os.mkdir( dportho )
            
        if os.path.isdir ( dres ) == False :
            os.mkdir( dres )
            
        # and create fasta files one per specie in one folder per group
        lsp = set( [ ] )
        fasta_files = { }
        
        for spprot in subfamily :
            sp,prot = spprot.split( ";" )
            lsp.add( sp )
            fasta = ddfasta[ int( sp ) ][ prot ]
            pathfasta = os.path.join( dportho, "sp" + str( sp ) + ".fasta" )
            ret = fasta_files.setdefault( pathfasta , [ ] ).append( ">" + prot + "\n" + fasta + "\n" )
            
        for fasta_file in fasta_files :
            sortie = file( fasta_file ,"w" )
            
            for line in fasta_files[ fasta_file ] :
                sortie.write( line )
                
            sortie.close( )

        if len(lsp) < 2 : # cannot run protein ortho on a single specie cluster
            continue   
        
        pathlist = os.path.join( params.portho_dir , "list_sub"+str(cnt)+".dat" )
        listsub.append( (pathlist, cnt) )
        sortie = file( pathlist ,"w")
        
        for sp in lsp :
            pathfasta = os.path.join( dportho, "sp" + str( sp ) + ".fasta" )
            sortie.write( pathfasta+"\n" )
            
        sortie.close( )
        ##########################################################################
        # run proteinortho on each of this sub fasta species    
        # run only blast
        # a lot of proteinortho session, may be put a limit on number of processus      
        portho_path_log = os.path.join( dres , "info_proteinortho_" + str( cnt ) + ".log" )
        portho_path_out = os.path.join( dres , "info_proteinortho_" + str( cnt ) + ".dat" )
        cmd = params.path_proteinortho + " -dir=" + dportho 
        cmd += " -log="+portho_path_log+" -o="+portho_path_out+" "
        # add custom parameters if exist
        for param in params.portho_params : 
            cmd += param+" "
            
        # add the fata list
        cmd += pathlist 
        
        subdir = os.path.join( params.portho_dir, "sub"+str(cnt) )
        result = pool_blastonly.apply_async( subprocess_threaded_blastonly , ( cmd, subdir ) )
        dicres[ cnt ] = ( result, portho_path_out )

    pool_blastonly.close( )
    pool_blastonly.join( )
    
    for cnt in dicres :
        result, path = dicres[cnt]
        if  result.ready( ) != True and result.successful() != True :
            msg = "\nError : Unable to run blast for file {} ".format( path )
            error_clean( msg, 1, os.path.join(p.workdir,".lock"), p.verbose, starting_time )
            raise ProteinorthoError( msg ) 
        err = result.get( )
        if err != "" :
            error_clean( err, 1, os.path.join(p.workdir,".lock"), p.verbose, starting_time )
            raise ProteinorthoError( err )
        
    return listsub

def compute_similarity( luniq_da, tmpdir, p, starting_time, path_compute_similarity, sub_file_size = 100000 ) :
    """Compute the similarity between unique domain arrangement
    
    Parameters
    ----------
    luniq_da : list 
        the list of unique domain arrangement
    tmpdir : string
        path to the tmpdir directory
    p : argparse.Namespace
        software input parameters from argparse module
    starting_time : int 
        the starting time of the software
    path_compute_similarity : string 
        path to compute_similarity software
    sub_file_size : int 
        number of pairwise similarity to compute per thread

    Returns
    -------
    GDA : graph 
        networkx undirected Graph object, nodes are domain arrangement, edges
        are set between nodes if similarity between domain arrangement is 
        superior to p.cutoff parameter
    da_similarity : dict 
        a dictionary of dictionary d[da1][da2] = sim
    
    """    
    cnt = 0 
    cnt_file = 0    
    if p.verbose :
        timestamp( "... froms scratch ", starting_time )
    da_similarity = { }
    # GDA is a graph of domain arrangement similarity, node = domain arrangement, edge = similarity
    GDA = nx.Graph( )
    lp = [ ]
    da_dir_tmp = os.path.join( tmpdir , "tmp" )
    if os.path.isdir ( da_dir_tmp ) == False:
        os.mkdir ( da_dir_tmp ) 

    da_dir_tmp_out = os.path.join( da_dir_tmp, "list_da_"+str(cnt_file)+".dat" )
    da_list_out = file( da_dir_tmp_out ,"w" )
    pool_sim = multiprocessing.Pool( processes=p.nb_job ) 
    dres = { }
    nb_uniq_da = len( luniq_da )
    # for domain arrangement compute the pairwise similarity 
    for i, str_dai in enumerate( luniq_da ):
        for j in range( i, nb_uniq_da ) :
            str_daj = luniq_da[j]
            # enougth pairwise similarity ?
            if cnt == sub_file_size :
                # close file
                da_list_out.close( )
                # output result file 
                pathout_sim = os.path.join( da_dir_tmp, "list_da_"+str(cnt_file)+".dat" )
                # TODO nice ionice is temporary, to remove
                #cmd = "nice -19 ionice -c3 "+path_compute_similarity
                cmd = path_compute_similarity
                cmd += " -i "+pathout_sim
                cmd += " -c "+str(p.cutoff)
                cmd += " -m "+p.matrix
                cmd += " -o %d -w F "%p.order
                # add to pool, run on previous file
                result = pool_sim.apply_async( subprocess_threaded_sim, (cmd,pathout_sim) )
                dres[cnt_file] = (result, pathout_sim )
                cnt = 0
                cnt_file += 1 
                # create new pairwise file
                da_dir_tmp_out = os.path.join( da_dir_tmp, "list_da_"+str(cnt_file)+".dat" )
                da_list_out = file( da_dir_tmp_out , "w" )
            # write pairwise similarity to compute
            da_list_out.write( str_dai+" "+str_daj + "\n")
            # initialise graph node ...
            # and dictionary keys
            ret = da_similarity.setdefault( str_dai, { } )
            ret = da_similarity.setdefault( str_daj, { } )
            cnt += 1 
    # don't forget to run the last one
    if cnt != sub_file_size :
        da_list_out.close( )
        pathout_sim = os.path.join( da_dir_tmp, "list_da_"+str(cnt_file)+".dat" )
        # TODO nice ionice is temporary, to remove
        #cmd = "nice -19 ionice -c3 "+path_compute_similarity
        cmd = path_compute_similarity
        cmd += " -i "+pathout_sim
        cmd += " -c "+str(p.cutoff)
        cmd += " -m "+p.matrix
        cmd += " -o %d -w F "%p.order
        # add to pool, run on previous file
        result = pool_sim.apply_async( subprocess_threaded_sim, (cmd,pathout_sim) )            
        dres[cnt_file] = (result, pathout_sim )
    # close and join pool
    pool_sim.close( )
    pool_sim.join( )
    # check each result and gather similarity of each pairwise DA compute
    path_lock = os.path.join( p.workdir, ".lock" )
    for cnt in dres :
        result, path_target = dres[ cnt ]
        # for some reason the thread is not ready after join/close
        if result.ready( ) != True :
            msg = "Error : Similarity for file {} is not finish\n".format(path_target)
            msg +="        multiprocessing Pool join and close function may be not working properly"
            error_clean( msg, 1, path_lock, p.verbose, starting_time )
            raise DASimilarityError( msg )
        # the thread is not successful at computing similarity
        if result.successful() != True :
            msg = "Error : Problem in running similarity for file {} ".format(path_target)
            msg += "        multiprocessing Pool join and close function may be not working properly"
            msg += "{}".format( result.get( ) )
            error_clean( msg, 1,path_lock, p.verbose, starting_time )
            raise DASimilarityError( msg )
        # get results
        out, err = result.get()
        if err != '' :
            msg = "Error : Problem in running  {} ".format(err)
            error_clean( msg, 1, path_lock, p.verbose, starting_time )
            raise DASimilarityError( msg )
        # parse the output result
        for line in out.split("\n") :
            if line == "" : continue 
            if line[0] == "\n" : continue 
            tmp = line.split( )
            try :
                sim = float( tmp[2] ) # change compute_similarity to get cosine or mwm
            except :
                msg = "Error : raising exeption when parsing compute_similarity results"
                error_clean( msg, 1, path_lock, p.verbose, starting_time )
                raise DASimilarityError( msg )
            # construct similarity graph and dictionary
            da_similarity[ tmp[0] ][ tmp[1] ] = sim
            GDA.add_edge( tmp[0], tmp[1], weight = sim )
    if p.verbose :
        timestamp( "done", starting_time )
    # get DA without edges
    missing = [ ]
    for da in luniq_da :
        if not da in GDA :
            missing.append( da )
    return GDA, da_similarity, missing

def cluster_domains( GDA, missing_da,  p, starting_time) :
    """Clusters the domain graph using DBSCAN algorithm and make a picture of
    the whole matrix map
    
    Parameters
    ----------
        GDA : Graph
            undirected graph of domain similarities
        missing_da : list
            DA not in Graph, no edge (no similarity) no any other DA
        p : argument parser object
            parameter object
        starting_time : int 
            program starting time
    Returns
    -------
        clusters : list
            a list of list containing the clusterised DA
    """
    # add missing da as a self cluster
    clusters = [ [da] for da in missing_da ]       
    
    if p.daonly :
        #if True :
        # to gain some memory space DBSCAN is only used on connected components
        clusters_comp = nx.connected_components( GDA )
        all_unclustered = [ ]
        for comp in clusters_comp :
            if len( comp ) > p.minpts :
                H = GDA.subgraph( comp )
                # networkx return an numpy.matrixlib.defmatrix.matrix 
                mat =  1.0 - np.array(nx.to_numpy_matrix( H, nodelist=comp ) ) 
                mat.flat[ :: mat.shape[0] + 1 ] = 0  # diag to 0
                # run OPTICS on distance matrix
                optics = Optics(  p.minpts, epsilon=p.epsilon )
                ordered, reachability, core_dist = optics.run(mat)
                labels = optics.cluster( p.epsilon_p ) 
                # run dbscan on distance matrix
                slabels = set( labels )
                for k in slabels :
                    ind = np.where( labels == k )[0]
                    if k == -1 :
                        for i in ind :
                            clusters.append( [comp[i]] )
                    else :
                        clusters.append(  [ comp[i] for i in ind ]  )
            else :
                # if the component is whith less memebers than the minpts cutoff
                # all the members of the same components are put in the same clusters
                clusters.append( comp )
    else :
        nodes = GDA.nodes( )
        # networkx return an numpy.matrixlib.defmatrix.matrix 
        # instead of an numpy.ndarray matrix, not really convenient ...
        bigmat = 1.0 - np.array( nx.to_numpy_matrix( GDA, nodelist=nodes ) )
        bigmat.flat[ :: bigmat.shape[0] +1 ] = 0 # diagonal to 0
        clusterid, error, nfound = Pycluster.kmedoids (bigmat, nclusters=p.kcluster, npass=10 )
        for l in np.unique( clusterid ) :
            tmp_clust =  [nodes[i] for i in range(clusterid.shape[0]) if clusterid[i] == l ] 
            clusters.append( tmp_clust )            

    return clusters
        