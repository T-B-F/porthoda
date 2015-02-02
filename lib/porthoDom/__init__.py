__version__ = '0.5'
__release__ = __version__ # + '-dev' # comment out '-dev' before a release

import traceback, inspect
import sys
import proteinorthoDom as porthoDom

__all__ = porthoDom.__all__


def porthoDom_main():
    porthoDom.porthoDom_main( )

if __name__ == "__main__" :
    """
    Catching malfunctionning behavior
    """
    try :
        ext = porthoDom_main( ) 
    except (LockError, ArgumentError, ProteinorthoError, DASimilarityError, ResultError, DependencyError, ExecutionError) as e :
        print "Software execution error : "
        print e.value
        traceback.print_exc()
        ext = 1 
    except SystemExit as e :
        ext = e.code
    except :
        print "Unexpected error:", sys.exc_info()[0]
        traceback.print_exc()
        ext = 1 
    finally :    
        sys.exit( ext )

    
