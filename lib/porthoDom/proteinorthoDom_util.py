﻿# -*- coding: utf-8 -*-
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

__author__ = "Tristan Bitard Feildel"
__email__ = "t.bitard.feildel@uni-muenster.de"
__institute__ = "Insitute for Evolution and Biodiversity"
__lab__ = "Evolutionary Bioinformatics"
__vesion__ = "1.0"

import ConfigParser, os, time, sys
from porthoDom.proteinorthoDom_err import ArgumentError, LockError
import shutil, multiprocessing
import networkx as nx

__all__ = [ "should_wait", "check_program", "timestamp", "error_clean", 
            "remove_dir", "last_dir" ]

# checking that the environment variable is present
if not os.environ.get('PORTHODOM'):
    print >>sys.stderr, ("Please set the PORTHODOM environ variable to the "
                         "location of PATH.ini")
    sys.exit(1)

config = ConfigParser.ConfigParser()
pathconfig=os.path.join(os.getenv("PORTHODOM"),"PATH.ini")
config.read(pathconfig)


def last_dir( path ) :
    """ Return the last directory of a path
    
    Parameter
    ---------
    path : string
        directory path
        
    Return
    ------
    tail : string
        last directory
        
    >>> p1 = "/toto/tatat/"
    >>> p2 = "toto/tatat/"
    >>> p3 = "/toto/tatat"
    >>> p4 = "toto/tatat"
    >>> p5 = "tatat"
    >>> p6 = ""
    >>> last_dir( p1 )
    'tatat'
    >>> last_dir( p2 )
    'tatat'
    >>> last_dir( p3 )
    'tatat'
    >>> last_dir( p4 )
    'tatat'    
    >>> last_dir( p5 )
    'tatat'
    >>> last_dir( p6 )
    ''
    """
    head, tail = os.path.split(path) 
    if tail == "" and head != "" :
            head,tail = os.path.split(head)
    elif tail == "" :
            tail = path 
    return tail

def remove_dir( path ) :
    """Check a path and remove a directory path
    Parameter
    ---------
    path : string
        the potential directory path
    """
    if os.path.isdir( path ) :
        shutil.rmtree( path )

def error_clean( msg, status, path_lock, verbose, starting_time ) :
    """Display an error message and exit after cleaning the lock
    
    Parameters
    ----------
    msg string
        the error message to display
    status int
        integer code for the program error 
    path_lock string
        path to the .lock file
    verbose bool
        display or not in verobse mode
    starting_time int
        starting time of the software (needed by timestamp function )
    """
    print >>sys.stderr, msg
    if os.path.isfile( path_lock ) :
        os.remove( path_lock )
    if verbose :
        timestamp( "Exiting", starting_time )

    #sys.exit( status )

################################################################################
# Utilitaries functions

class storage :
    """
    Empty class, just use to store parameters and give them to a function
    (Less parameters in the functions header, even if it's not 'clean')
    """
    pass

def should_wait( lp, time2sleep, size ) :
    """ A waiting function, to manage and to wait for subprocess worker to terminate
    
    Parameters 
    ----------
        lp : list 
            the list of process
        time2sleep : int 
            number of second to wait if maximal size is reached
        size : int 
            maximal number of worker
            
    Returns
    -------
        lp : list 
            the list of process without the one terminated
        processed : list
            the list of processed finished
        lout : list 
            the list of output message form terminated process
        lerr : list 
            the list of error message form terminated process
    """
    lout = [ ]
    lerr = [ ]
    processed = [ ]
    while len(lp) >= size :
        to_remove = [ ]
        time.sleep( time2sleep )
        for process in lp:
             if process.poll( ) != None : # finish ?
                 out,err = process.communicate( )   
                 lout.append( out )
                 lerr.append( err )
                 processed.append( process )
                 to_remove.append( process ) 
        for process in to_remove :
            lp.remove( process )
    return lp, processed, lout, lerr

def check_program(program):
    """ Check in the $PATH if we found the program
    from http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    
    Parameter
    ---------
        program : string
            the program name to look for
    Return
    ------
        path : string 
            the path to the program or None if the program is not found
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    def ext_candidates(fpath):
        yield fpath
        for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
            yield fpath + ext

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            for candidate in ext_candidates(exe_file):
                if is_exe(candidate):
                    return candidate
    return None
    
def check_program_config(sections=[]):
    """ check the program of the configparser
    
    Parameter
    ---------
    sections: list
        list of sections to check
    Return
    ------
    allhere : bool
       True if all programs are found, else False
    """
    global config
    allhere = True
    if not sections:
        sections = config.sections()
    for section in sections:
        for key, program in config.items(section):
            path = check_program(program)
            print section, key, program, path
            if path == None:
                print >>sys.stderr, ("Error : unable to find program path {} "
                                     "from config file ".format(program))
                allhere=False
    return allhere
        
def timestamp(s,start_time):
    """ Display the ealpsed time since the start of the program with a specific
    sentence
    
    Parameter
    ---------
    s : string 
        sentence to display
    start_time : int 
        the starting time of the program
    """
    t = time.gmtime(time.time() - start_time)
    print s.ljust(44), time.strftime('%H:%M:%S', t)  


def check_parameters(params, path_lock, starting_time):
    """ check the compatibility between the different parameters, raise an 
    ArgumentError if an incompatibility occurs
    
    Parameters
    ----------
    params : object
        argparse object argument parser
        
    Return :
    bothpfamonlyanddone: bool
        True if both pfamonly and pfamdone parameters are set to False
    bothblastonlyanddone: bool
        True if both blastonly and blastdone parameters are set to False
    """
    # number of cpu
    nbcpu = multiprocessing.cpu_count()
    if params.nb_job > nbcpu:
        msg = "number of cpu requested {} ".format(params.nb_job)
        msg += "superior to the number of cpu present {}".format(nbcpu)
        if params.force == False:
            msg = "\nError: "+ msg
            error_clean(msg, 1, path_lock, params.verbose, starting_time)
            raise ArgumentError(msg)
        else:
            print >>sys.stderr, "\nWarning: "+msg
            
    if params.pfamonly and params.pfamdone:
        msg = "\nError: Cannot use both pfamonly and pfamdone at the same time"
        msg += "\n"
        error_clean(msg, 1, path_lock, p.verbose, starting_time)
        raise ArgumentError(msg)
    
    if params.pfamonly == False and params.pfamdone == False:
        bothpfamonlyanddone = True
    else:
        bothpfamonlyanddone = False

    # parameter check
    if params.blastdone and params.blastonly:
        msg = "\nError: Cannot use both blastonly and blastdone at the same "
        msg += "time\n"
        error_clean(msg, 1, path_lock, params.verbose, starting_time)
        raise ArgumentError(msg)

    if params.blastonly == False and params.blastdone == False:
        bothblastonlyanddone = True
    else:
        bothblastonlyanddone = False
        
    if params.daonly and (params.blastonly or params.blastdone):
        msg = "\nError: Cannot use domain arrangement similarity and sequence"
        msg += "similarity clustering at the same time\n"
        error_clean(msg, 1, path_lock, params.verbose, starting_time)
        raise ArgumentError(msg)
    
    return bothpfamonlyanddone, bothblastonlyanddone

def prepare_lock(params, path_lock, path_lock_param):
    """  check lock files and parameters for lock, raise an LockError if any
    problem"
    
    Parameters
    ----------
    params : object
        the argparse argument parser object
    path_lock: string
        the path to the current lock
    pathlock_param: string
        the path to the file storing the parameters of the current working 
        directory
    """
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
        old_cutoff = params.cutoff
        olf_order = params.order
        with open(path_lock_param) as f:
            for i, line in enumerate(f):
                if i == 0:
                    old_cutoff = float(line.strip())
                elif i == 1:
                    old_order = int(line.strip())
            if old_cutoff != params.cutoff or old_order != params.order:
                msg = u"""
Error: Previous result of this directory doesn't match
       new input parameters for domain similarity
       Please use the same parameter or use a different
       working directory.
       Old parameters: cutoff {}, order {}.
       New parameters: cutoff {}, order {}.
       Proteinortho parameters can't be changed without any problems,
       specify new parameters with --portho parameters.""".format(old_cutoff, 
       old_order, params.cutoff, params.order)
                raise LockError(msg)
            else:
                if params.verbose:
                    print ("Same previous parameters detected for domain "
                        "similarity, continue")    
