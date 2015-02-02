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
""" Error classes for different context
"""

import os, sys

__author__ = "Tristan Bitard Feildel"
__email__ = "t.bitard.feildel@uni-muenster.de"
__institute__ = "Insitute for Evolution and Biodiversity"
__lab__ = "Evolutionary Bioinformatics"
__vesion__ = "1.0"

__all__ = ["MsgError","LockError","ArgumentError","ProteinorthoError",
"DASimilarityError","ResultError", "DependencyError", "ExecutionError" ]
    
    
class MsgError( Exception ) :     
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return str(self.value)

class LockError( MsgError ) :
    def __init__( self, value) :
        MsgError.__init__( self, value )
       
class ArgumentError( MsgError ) :
    def __init__( self, value) :
        MsgError.__init__( self, value )

class ProteinorthoError( MsgError ) :
    def __init__( self, value) :
        MsgError.__init__( self, value )

class DASimilarityError( MsgError ) :
    def __init__( self, value) :
        MsgError.__init__( self, value )
        
class ResultError( MsgError ) :
    def __init__( self, value) :
        MsgError.__init__( self, value )

class DependencyError( MsgError ) : 
    def __init__( self, value) :
        MsgError.__init__( self, value )

class ExecutionError( MsgError ) :
    def __init__( self, value) :
        MsgError.__init__( self, value )
      