# -*- coding: utf-8 -*-
#!/usr/bin/env python
# utilitaries functions for proteinorthoDA.py
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

"""
 proteinorthoDA is a wrapper around proteinortho.pl script using a 
 pre-clustering step based on domain arrangement similarity to speed up the 
 blast all against all step used by proteinorhto.
"""

############################################################################
# suggested tree structure of the program directories output
#
# workdir /
#     |
#     | .lock
#     | ----- da /
#     | ----- tmpdir (exp 1)/
#                | .lock_param
#                | ----- fasta / 
#                | ----- portho / 
#                | ----- result / 
#
# or 
# workdir /
#     |
#     | .lock
#     | ----- da /
# tmpdir (exp 1)/
#     | .lock_param
#     | ----- fasta / 
#     | ----- portho / 
#     | ----- result / 

__author__ = "Tristan Bitard Feildel"
__email__ = "t.bitard.feildel@uni-muenster.de"
__institute__ = "Insitute for Evolution and Biodiversity"
__lab__ = "Evolutionary Bioinformatics"
__vesion__ = "0.5"

import os, sys, subprocess, multiprocessing, pickle
import time, tarfile, shutil, shlex, argparse


from porthoDA.proteinorthoDA_err import LockError, ArgumentError, ResultError
from porthoDA.proteinorthoDA_err import ProteinorthoError, DASimilarityError
from porthoDA.proteinorthoDA_err import ExecutionError, DependencyError 
from porthoDA.proteinorthoDA_util import should_wait, check_program, timestamp
from porthoDA.proteinorthoDA_util import error_clean, remove_dir, storage
from porthoDA.proteinorthoDA_util import check_program_config
from porthoDA.proteinorthoDA_worker import subprocess_threaded_blastdone
from porthoDA.proteinorthoDA_IO import read_multifasta, read_porthoparams
from porthoDA.proteinorthoDA_IO import write_results_daclusters, extractPFidDA
from porthoDA.proteinorthoDA_IO import write_results_porthodaclusters
from porthoDA.proteinorthoDA_IO import read_dadone
from porthoDA.proteinorthoDA_algo import run_proteinortho_blast
from porthoDA.proteinorthoDA_algo import cluster_domains, compute_similarity

__all__ = ['porthoDA_main']

def process_arguments()
    """ process arguments with the argument parser module
    
    Return
    ------
    params: object
        argparse object storing arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", action="store", dest="workdir",  
        help="working directory")
    parser.add_argument("-t", action="store", dest="tmpdir", default=None,
        help="temporary directory for proteinortho intermediate results")
    parser.add_argument("-i", action="store", dest="listproteomes",  
        help="files with proteomes list")
    parser.add_argument("-c", action="store", dest="cutoff", type=float,
        help="cutoff similarity, only score higher than cutoff are conserved")
    parser.add_argument("-m", action="store", dest="matrix", 
        help="path to matrix file (must correspond to the pfam database)")
    parser.add_argument("-p", action="store", dest="pfamdb", 
        help="path to pfam database (must correspond to the matrix)")
    parser.add_argument("-O", action="store", dest="order", default=1, type=int, 
       help="control the order of magnitude used by compute_similarity [def:1]")
    parser.add_argument("-o", action="store", dest="output", help="output file")
    parser.add_argument("-n", action="store", dest="nb_job", 
        help="number of job running at the same time", type=int, default=50)
    parser.add_argument("-k", action="store", dest="kcluster", type=int , 
        help="number of cluster")
    parser.add_argument("--maskrepeat", action="store_true", dest="maskrepeat",
        help="mask repeat in proteins", default=False) 
    parser.add_argument("-f", "--force", action="store_true", dest="force", 
        help="force parameters", default=False)
    parser.add_argument("--daonly", action="store_true", dest="daonly", 
        help="cluster only on domain arrangement similarity", default=False)
    parser.add_argument("--dadone", action="store_true", dest="dadone", 
        help="choose to not run domain similarity", default=False)
    parser.add_argument("--pfamonly", action="store_true", dest="pfamonly", 
        help="choose to only run pfamscan", default=False)
    parser.add_argument("--pfamdone", action="store_true", dest="pfamdone", 
        help="choose to not run pfamscan", default=False)
    parser.add_argument("--blastonly", action="store_true", dest="blastonly", 
        help="choose to only run blast", default=False)
    parser.add_argument("--blastdone", action="store_true", dest="blastdone", 
        help="choose to not run blast", default=False)        
    parser.add_argument("--portho", action="store", dest="porthoparams", 
        help="run proteinortho with specific parameters", default=None)
    parser.add_argument("--optics-eps", action="store", dest="epsilon", 
        help="The upper bound epsilon parameters for OPTICS algorithms", 
        type=float, default=0.5)
    parser.add_argument("--optics-epsprime", action="store", dest="epsilon_p", 
        help="The second epsilon parameters for OPTICS algorithms", type=float, 
        default=0.2)
    parser.add_argument("--optics-minpts", action="store", dest="minpts", 
        help="The minpoints parameters for OPTICS algorithms", type=float, 
        default=5)
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose", 
        help="verbose information", default=False)
    parser.add_argument("--clean-lock", action="store_true", dest="cleanlock", 
        help="clean the lock files", default=False)
    parser.add_argument("--clean-proteinortho", action="store_true", dest="cleanproteinortho", 
        help="clean the proteinortho directories (ie -t dir)", default=False)
    parser.add_argument("--clean-pfamscan", action="store_true", dest="cleanpfamscan", 
        help="clean the pfamscan directory", default=False)
    parser.add_argument("--clean-all", action="store_true", dest="cleanall", 
        help="clean the working directory", default=False)

    params = parser.parse_args()
    return params
    
# main fonction
def porthoDA_main():
    """ proteinorthoDA is a python package and software build around the 
    orthologous protein detector proteinortho.
    proteinorthoDA use domain annotation to restrict the orthologous proteins
    space search of proteinortho by only looking at protein sharing a similar
    domain arrangement as potential orthologous proteins.
    """
    
    ############################################################################
    # Check if lock file exist 
    
    if p.tmpdir == None:
        name = os.path.join(p.workdir, "exp_{}_{}_{}_{}_{}".\
            format(p.cutoff, p.order, p.epsilon, p.minpts, p.maskrepeat))
        if p.verbose: 
            print "Path of the temporary directory for experiment:\n\t"+name
        p.tmpdir = name 

    path_lock = os.path.join(p.workdir, ".lock")
    path_lock_param = os.path.join(p.tmpdir, ".lock_param")

    ############################################################################
    # clean all ? 
    
    if p.cleanall:
        if p.verbose: 
            print "Removing the working directory"
        remove_dir(p.workdir)
        if p.verbose: 
            print "done"
        return 0

    ############################################################################
    # clean pfamscan ? 
    
    if p.cleanpfamscan:
        if p.verbose: 
            print "Removing the pfamscan directory"
        da_dir = os.path.join(p.workdir, "da")
        remove_dir(da_dir)
        if p.verbose: 
            print "done"
        return 0            

    ############################################################################
    # clean proteinortho ? 
    
    if p.cleanproteinortho:
        if p.verbose: 
            print "Removing the proteinortho directories \
                                          (<tmpdir>/result and <tmpdir>/portho)"
        remove_dir(os.path.join(p.tmpdir, "result"))
        remove_dir(os.path.join(p.tmpdir, "portho"))
        if p.verbose: 
            print "done"
        return 0  
        

    ############################################################################
    # clean lock ? 
    
    if p.cleanlock:
        if p.verbose: 
            print "Removing lock files ..."
        if os.path.isfile(path_lock) == True:
            os.remove(path_lock)
        if os.path.isfile(path_lock_param) == True:
            os.remove(path_lock_param)
        if p.verbose: 
            print "done"
        return 0


    ############################################################################
    # creating lock
        
    if os.path.isfile(path_lock):
        msg = u"Error: This directorty is currently used or the job using "
        msg += u"this directory didn't reach the end of the script.\n"
        msg += u"        Please choose an other directory or delete the .lock* "
        msg += u"files inside the directory.\n"
        msg += u"        (you can use --clean-lock option for this operation)\n"
        msg += u"        If previous analysis have been done in this directory,"
        msg += u"they would disturbe the new one.\n"
        msg += u"        To clean proteinorhto result use --clean-proteinortho "
        msg += u"option.\n"
        msg += u"        To clean pfamscan result use --clean-pfamscan option."
        msg += u"\n        To clean everything use --clean-all.\n"
        raise LockError(msg)

    if os.path.isfile(path_lock_param): 
        with open(path_lock_param) as inf:
            old_cutoff = p.cutoff
            old_order = p.order
            for i, line in enumerate(inf):
                if i == 0:
                    old_cutoff = float(line.strip())
                elif i == 1:
                    old_order = int(line.strip())
        if old_cutoff != p.cutoff or old_order != p.order:
            msg = u"Error: Previous result of this directory doesn't match\n"
            msg += u"        new input parameters for domain similarity\n"
            msg += u"        Please use the same parameter or use a different "
            msg += u"working directory.\n"
            msg += u"        Old parameters: cutoff {}, order {} .\n"\
                .format(old_cutoff, old_order)
            msg += u"        New parameters: cutoff {}, order {} .\n"\
                .format(p.cutoff, p.order)
            msg += u"        Proteinortho parameters can't be changed without "
            msg += u"any problems, \n"
            msg += u"        specify new parameters with --portho parameters.\n"
            raise LockError(msg)
        else:
            if p.verbose:
                print ("Same previous parameters detected for domain "
                       "similarity, continue")

    # start here, after error message
    starting_time = time.time()
    if p.verbose: # only if verbose mose activated
        timestamp("Starting ", starting_time)
        

    ###########################################################################
    # programm check (pfamscan and proteinortho)        

    # number of cpu
    nbcpu = multiprocessing.cpu_count()
    if p.nb_job > nbcpu:
        msg = "number of cpu requested {} ".format(p.nb_job)
        msg += "superior to the number of cpu present {}".format(nbcpu)
        if p.force == False:
            msg = "\nError: "+ msg
            error_clean(msg, 1, path_lock, p.verbose, starting_time)
            raise ArgumentError(msg)
        else:
            print >>sys.stderr, "\nWarning: "+msg
            
    # check that program in PATH.init (read by config in proteinorthoDA_util.py
    # are all here
    if not check_program_config():
        msg = ("\nError: Unable to find all programs, please check your "
        "PATH.init")
        error_clean(msg, 1, path_lock, p.verbose, starting_time)
        raise DependencyError(msg)
    path_pfam_scan = config.get("annotation", "pfamscan")
    check_proteinortho = config.get("similarity", "dasim")
    
    if p.verbose:
        print path_pfam_scan.ljust(44), " ... found"

    # proteinortho.pl        
    check_proteinortho = check_program("proteinortho4.pl")
    if check_proteinortho == None:
        msg = "\nError: Unable to find proteinortho4.pl script"
        error_clean(msg, 1, path_lock, p.verbose, starting_time)
        raise DependencyError(msg)
    else:
        path_proteinortho = check_proteinortho
    if p.verbose:
        print path_proteinortho.ljust(44), " ... found"
    
    # check compute_similarity 
    check_compute_similarity = check_program("compute_similarity")
    if check_proteinortho == None:
        msg = "\nError: Unable to find compute_similarity software"
        error_clean(msg, 1, path_lock, p.verbose, starting_time)
        raise DependencyError(msg)
    else:
        path_compute_similarity = check_compute_similarity
    if p.verbose:
        print path_compute_similarity.ljust(44), " ... found"
        
    # check if user want some specific parameters for proteinortho
    if p.porthoparams != None:
        portho_params = read_porthoparams(p.porthoparams, path_proteinortho)
    else:
        portho_params = []

    ###########################################################################
    # parameter check
    if p.pfamonly and p.pfamdone:
        msg = "\nError: Cannot use both pfamonly and pfamdone at the same time"
        msg += "\n"
        error_clean(msg, 1, path_lock, p.verbose, starting_time)
        raise ArgumentError(msg)
    
    if p.pfamonly == False and p.pfamdone == False:
        bothpfamonlyanddone = True
    else:
        bothpfamonlyanddone = False

    # parameter check
    if p.blastdone and p.blastonly:
        msg = "\nError: Cannot use both blastonly and blastdone at the same "
        msg += "time\n"
        error_clean(msg, 1, path_lock, p.verbose, starting_time)
        raise ArgumentError(msg)

    if p.blastonly == False and p.blastdone == False:
        bothblastonlyanddone = True
    else:
        bothblastonlyanddone = False
        
    if p.daonly and (p.blastonly or p.blastdone):
        msg = "\nError: Cannot use domain arrangement similarity and sequence"
        msg += "similarity clustering at the same time\n"
        error_clean(msg, 1, path_lock, p.verbose, starting_time)
        raise ArgumentError(msg)

    ###########################################################################
    # create or check if working directory is empty 
    if os.path.isdir(p.workdir) == False:
        os.makedirs(p.workdir)
        
    if os.path.isdir(p.tmpdir) == False:
        os.makedirs(p.tmpdir)
        
    # create sub directories
    #fasta_dir = os.path.join(p.tmpdir, "fasta")
    da_dir = os.path.join(p.workdir, "da")
    portho_dir = os.path.join(p.tmpdir, "portho")
    res_dir = os.path.join(p.tmpdir, "result")
    #if os.path.isdir (fasta_dir) == False:
    #    os.makedirs(fasta_dir)
    if os.path.isdir(da_dir) == False:
        os.makedirs(da_dir)
    if os.path.isdir(res_dir) == False:
        os.makedirs(res_dir)
    if os.path.isdir(portho_dir) == False:
        os.makedirs(portho_dir)

    ###########################################################################
    # Create lock file and lock parameter
    lock = file(path_lock, "w")
    lock.close()
    lock = file(path_lock_param, "w")
    lock.write("{}\n".format(p.cutoff))
    lock.write("{}\n".format(p.order))
    lock.close()        
        
    ###########################################################################
    # run pfams_can.pl
    
    # read file with proteomes.fasta path
    path_proteomes = {}
    name_proteomes = {}
    ddfasta = {}
    cnt = 0    
    with open(p.listproteomes) as inf:
        if p.verbose:
            timestamp("Reading fasta files ... ", starting_time)
        for line in inf:
            path_proteomes[cnt] = line.strip()
            _, name_proteomes[cnt] = os.path.split(line.strip())
            dfasta = read_multifasta(line.strip())
            ddfasta[cnt] = dfasta
            cnt += 1 
    if p.verbose:
        print "done"

    # run pfamcan in parallel
    if bothpfamonlyanddone or p.pfamonly:
        if p.verbose:
            timestamp("Running pfam scan analysis ...", starting_time)
        lp = []
        for k in path_proteomes:
            pfam_scan_out = os.path.join(da_dir, str(k)+".pfamscan")
            cmd = "{} -fasta {} -dir {}, -outfile {} ".format(
                path_pfam_scan, path_proteomes[k], p.pfamdb, pfam_scan_out)
            cmd = shlex.split(cmd)
            pid = subprocess.Popen(cmd)
            lp.append(pid)
            lp, processed, lout, lerr = should_wait(lp, 3, p.nb_job) 
        # wait until all proteomes are finish
        for pid in lp:
            pid.wait()
        if p.verbose:
            timestamp("done", starting_time)
    # only pfamscan running ?
    if p.pfamonly:
        os.remove(path_lock)
        sys.exit(0)

    ###########################################################################
    # read results and get unique da

    if p.pfamdone: # check output directory of pfam 
        l = os.listdir(da_dir) 
        if l == []:
            msg = "\nError: pfamscan output directory is empty\n"
            msg += "        please check the working directory for pfamscan "
            msg += "output file or run the program with --pfamonly or without "
            msg += "--pfamdone options\n"
            error_clean(msg, 1, path_lock, p.verbose, starting_time)       
            raise ExecutionError(msg)

    uniq_da = {}
    sp_without_DA = {}
    for k in path_proteomes:
        sp_without_DA[k] = []
        pfam_scan_out = os.path.join(da_dir, str(k)+".pfamscan")
        d = extractPFidDA(pfam_scan_out, mask_repeat=p.maskrepeat)
        fastas = ddfasta[k]
        for prot in fastas:
            # store proteins without domain arrangement
            if d.has_key(prot) == False:
                sp_without_DA[k].append(prot)
            # store proteins with domain arrangement
            else:
                da = ";".join(d[prot])
                # one DA associated with a list of sp;prot
                ret = uniq_da.setdefault(da, []).append(str(k)+";"+prot) 

    ###########################################################################
    # compute DA similarity
    
    # compute similarity
    luniq_da = uniq_da.keys()
    nb_uniq_da = len(luniq_da)
    if p.verbose:
        print "Number of unique DA {} (with masked repeat: {}) "\
        .format(nb_uniq_da, p.maskrepeat)

    if not p.daonly:
        if p.kcluster == None:
            if p.verbose:
                print "Automaticaly adjust the number of cluster"
            # how many proteins per clusters ? -> 20
            p.kcluster = int(nb_uniq_da / 20) 
            if nb_uniq_da > p.kcluster * 20:
                p.kcluster += 1 # round to upper boundarie
            if p.kcluster < 20:
                p.kcluster = 20 

        if p.verbose:
            print ("Number of cluster used to diminue the sequence "
                   "space search: {}".format(p.kcluster))
        
    if p.daonly:
        if not p.dadone:
            try:
                GDA, da_similarity, missing_da = compute_similarity(
                luniq_da, p.tmpdir, p, starting_time, path_compute_similarity)
            except DASimilarityError:
                raise DASimilarityError
            except:
                raise
            
            sortie = file(os.path.join(p.tmpdir, "da_dict_sim.dat"), "w")
            pickle.dump(da_similarity, sortie) 
            sortie.close()    
        else:
            GDA, da_similarity, missing_da = read_dadone(
                            os.path.join(p.tmpdir, "da_dict_sim.dat"), luniq_da)
        
        ########################################################################
        # gather DA, clustering method ? connected component after cutoff
        if p.verbose:
            timestamp("Extract clusters", starting_time)
        
        clusters = cluster_domains(GDA, missing_da, p, starting_time)
        
        del GDA
        
        if p.verbose:
            timestamp("done", starting_time)
            timestamp("Writing results", starting_time)            
        
        write_results_daclusters(clusters, uniq_da, name_proteomes, p)
        
        if p.verbose:
            timestamp("done", starting_time)
        if p.verbose:
            timestamp("Finished", starting_time)
        os.remove(path_lock)
        # End of clustering DA only
        sys.exit(0)
    
    if bothblastonlyanddone or p.blastonly:
        if p.verbose:
            print "Computing similarity ... ",
            
        if p.dadone:
            if p.verbose:
                timestamp("reading dictionary", starting_time)
            GDA, da_similarity, missing_da = read_dadone(
                            os.path.join(p.tmpdir, "da_dict_sim.dat"), luniq_da)
                
        else:
            # output temporary sub file of domain pair similarity to be computed
            # each file correspond to a subprocess job
            try:
                GDA, da_similarity, missing_da = compute_similarity(
                luniq_da, p.tmpdir, p, starting_time, path_compute_similarity)
            except DASimilarityError:
                raise DASimilarityError
            except:
                raise
                        
            sortie = file(os.path.join(p.tmpdir, "da_dict_sim.dat"), "w")
            pickle.dump(da_similarity, sortie) 
            sortie.close() 
        ########################################################################
        # gather DA, clustering method ? connected component after cutoff

        clusters = cluster_domains(GDA, missing_da, p, starting_time)   
    
        del GDA
    
    ###########################################################################
    # write sub fasta files corresponding to gathered DA / divid per species
        if p.verbose:
            print "number of clusters", len(clusters)    
    
    params = storage()
    params.verbose = p.verbose
    params.portho_dir = portho_dir
    params.res_dir = res_dir
    params.nb_job = p.nb_job
    params.portho_params = portho_params
    params.path_proteinortho = path_proteinortho
    
    # blastonly and blastdone
    if bothblastonlyanddone:        
        if p.verbose:
            timestamp("run proteinorho", starting_time) 
            
        try:
            listsub = run_proteinortho_blast(clusters, uniq_da, ddfasta, params)     
        except ProteinorthoError as error:
            raise error
        except:
            raise
        
        if p.verbose:
            timestamp("done", starting_time)
            
    # blastonly
    elif p.blastonly:
        if p.verbose:
            timestamp("run proteinorho mode blastonly", starting_time)
            
        params.portho_params.appen("-blastonly")
        try:
            listsub = run_proteinortho_blast(clusters, uniq_da, ddfasta, params)
        except ProteinorthoError as error:
            raise error
        except:
            raise
        
        if p.verbose:
            timestamp("done", starting_time)
        # end of blatonly
        sys.exit(0)
        
    # blastdone
    elif p.blastdone:
        listsub = []
        
        if p.verbose:
            timestamp("run proteinortho mode blastdone ", starting_time) 
            
        tmpl = os.listdir(portho_dir)
        
        for tmpf in tmpl:
            if tmpf[-4:] == ".dat" and tmpf[:8] == "list_sub":
                cnt = int(tmpf[8:-4])
                path_tmpf = os.path.join(portho_dir, tmpf)
                listsub.append((path_tmpf, cnt))
                
        # run result analysis
        dpid = {}
        dresults = {}
        lp = []
        pool_blastdone = multiprocessing.Pool(processes=p.nb_job)
        
        for pathlist, cnt in listsub:    
            dportho = os.path.join(portho_dir, "sub" + str(cnt))
            dportho_gz = os.path.join(portho_dir, "sub" + str(cnt) + ".tar.gz")
            dres = os.path.join(res_dir, "sub" + str(cnt))
            dres_gz = os.path.join(res_dir, "sub" + str(cnt) + ".tar.gz")
            
            if os.path.isdir(dres) == False:
                os.mkdir(dres)
                
            if os.path.isfile(dres_gz):
                os.remove(dres_gz)
                
            if os.path.isfile(dportho_gz):
                tar = tarfile.open(dportho_gz, "r:gz")
                tar.extractall(path=portho_dir)
                tar.close()
            
            portho_path_log = os.path.join(dres, 
                                    "info_proteinortho_" + str(cnt) + ".log")
            portho_path_out = os.path.join(dres, 
                                    "info_proteinortho_" + str(cnt) + ".dat")
            cmd = path_proteinortho + " -dir=" + dportho + " -blastdone "
            cmd += "-log="+portho_path_log + " -o="+portho_path_out+" " #-a=2 "
            
            # add custom parameters if exist 
            for param in portho_params: 
                cmd += param + " "
            
            # add fasta file
            cmd += pathlist
            result = pool_blastdone.apply_async(
                            subprocess_threaded_blastdone, (cmd, dportho, dres)) 
            dresults[cnt] = (result, dportho)
            
        pool_blastdone.close()
        pool_blastdone.join()
        
        for cnt in dresults:
            result, path = dresults[cnt]
            
            if result.ready() != True or result.successful() != True:
                msg = "\nError: Unable to run blastdone on file {}".format(path)
                error_clean(msg, 1, path_lock, p.verbose, starting_time)
                raise ProteinorthoError(err)
            err = result.get()
            if err != "":
                error_clean(err, 1, path_lock, p.verbose, starting_time)
                raise ProteinorthoError(err)

        if p.verbose:
            timestamp("done", starting_time)
            
    ###########################################################################
    # gathered proteinortho results
    if p.verbose:
        timestamp("Gather proteinortho results", starting_time)
        
    # read results
    results = set([])
    dresults = {}
    portho_info = {}
    lsp = sorted(name_proteomes.keys())
    offset = 3 # proteinortho output starting position of informations
    fam = 0
    
    for pathlist, cnt in listsub: 
        dres = os.path.join(res_dir, "sub"+str(cnt))
        dres_gz = os.path.join(res_dir, "sub"+str(cnt)+".tar.gz")
        
        if os.path.isfile(dres_gz):
            tar = tarfile.open(dres_gz, "r:gz")
            tar.extractall(path=res_dir)
            tar.close()
            os.remove(dres_gz)
            
        elif os.path.isdir(dres) == False:
            msg = "\nError: Unable to find {} or {} ".format(dres_gz, dres)
            error_clean(msg, 1, path_lock, p.verbose, starting_time)
            raise ResultError(msg)
            
        path_res = os.path.join(dres, "info_proteinortho_"+str(cnt)+".dat")
        
        with open(path_res) as inf:
            line = inf.readline()
            # the format is sp<int>.fasta
            species = [int(sp[2:].split(".")[0]) for sp in line.split()[offset:]] 
            line = inf.readline() # second header
            # store every clusters
            for line in inf: # two header lines
                if line[0] == "#": 
                    footer = line 
                    continue
                
                else:
                    tmp = line.split()
                    nbspecies = tmp[0]
                    nbproteins = tmp[1]
                    nbconn = tmp[2]
                    prots = tmp[3:]
                    tmp = tuple()
                    
                    for i, prot in enumerate(prots):
                        if prot != "*":
                            tmp += ((species[i], prot),)
                            portho_info.setdefault(species[i], {})
                            
                            for p_ in prot.split(","):
                                portho_path_out = os.path.join(
                                    dres, "info_proteinortho_" + str(cnt) + ".dat")
                                if portho_info[species[i]].has_key(p_):
                                    print >>sys.stderr, "Possible multiple key entry for protein "+p_+" path: "
                                    print >>sys.stderr, "    "+portho_info[species[i]][p_]
                                    print >>sys.stderr, "    "+portho_path_out
                                    
                                portho_info[species[i]][p_] = portho_path_out
                                
                    dresults[fam] = (nbspecies, nbproteins, nbconn, tmp)
                    fam += 1
        
        tar = tarfile.open(dres_gz, "w:gz")
        tar.add(dres, arcname="sub"+str(cnt))
        tar.close()
        shutil.rmtree(dres)
        
    if p.verbose:
        timestamp("done", starting_time)
        timestamp("Write results", starting_time)
    write_results_porthodaclusters(name_proteomes, portho_info, dresults, 
                                   footer, p)

    ###########################################################################
    if p.verbose:
        timestamp("Exiting properly", starting_time)
    
    os.remove(path_lock)
    
    # exit successfuly 
    return 0
    
