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
""" Module to manage multiple worker threads for blast and similarity 
computation
"""

import os, sys
import subprocess
import shlex
import tarfile
import shutil
from proteinorthoDA_util import last_dir

__author__ = "Tristan Bitard Feildel"
__email__ = "t.bitard.feildel@uni-muenster.de"
__institute__ = "Insitute for Evolution and Biodiversity"
__lab__ = "Evolutionary Bioinformatics"
__vesion__ = "1.0"

__all__ = [ "subprocess_threaded_sim", "subprocess_threaded_blastonly", "subprocess_threaded_blastdone" ]
    
def subprocess_threaded_sim( command, path ) :
    """Function to run the subprocess similarity computataion
    
    Parameters:
    -----------
        command : string 
            command line argument to run
        path ; string 
            path to the file to remove after computataion
    
    Returns
    -------
        out : string 
            output message
        err : string 
            error message
    """
    if isinstance( command, basestring ) :
        command = shlex.split( command )
    cmd = command     
    try :
        out  = subprocess.check_output( cmd, universal_newlines=True )
        err = ""
    except subprocess.CalledProcessError, e :
        out = e.output
        err = "Error : running similarity : %s "%" ".join(e.cmd)
    os.remove( path )
    return out, err
    
def subprocess_threaded_blastonly( command, path ) :
    """Function to run the subprocess similarity computataion
    
    Parameters 
    ----------
        command : string 
            command line argument to run
        path : string 
            path to the file to remove after computataion
        
    Returns
    -------
        out : string 
            output message
        err : 
            string error message
    """    
    if isinstance( command, basestring ) :
        command = shlex.split( command )
    cmd = command 
    with open( os.devnull, 'w') as tempf:    
        try :
            ret  = subprocess.check_call( cmd, universal_newlines=True , stdout=tempf, stderr=tempf )
            err = ""
        except subprocess.CalledProcessError, e :
            err = "Error : running blast computation on\n\t{} ".format(" ".join(e.cmd))
    tar = tarfile.open(path+".tar.gz", "w:gz")
    tar.add( path, arcname=last_dir(path) )
    tar.close()
    shutil.rmtree( path )
    return err    
    
def subprocess_threaded_blastdone( command, dwork, dres ) :
    """ Function to run the subprocess similarity computataion
    
    Parameters
    ----------
        command : string 
            command line argument to run
        dword : string 
            path to the working directory 
        dres : string 
            path to the result directory
            
    Returns
    -------
        out : string 
            output message
        err : 
            string error message
    """    
    if isinstance( command, basestring ) :
        command = shlex.split( command )
    cmd = command     
    with open( os.devnull, 'w') as tempf:    
        try :
            ret  = subprocess.check_call( cmd, stdout=tempf, stderr=tempf,  universal_newlines=True )
            err = ""
        except subprocess.CalledProcessError, e :
            err = "Error : running blast analysis on\n\t{} ".format(" ".join(e.cmd))
    shutil.rmtree( dwork )
    tar = tarfile.open(dres+".tar.gz", "w:gz")
    tar.add( dres, arcname=last_dir(dres) )
    tar.close()
    shutil.rmtree( dres )
    return err   
    