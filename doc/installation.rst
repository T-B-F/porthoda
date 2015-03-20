.. _installation:

Installation
============

Requirement
-----------

:program:`porthoDom` is using several tools that need to be installed 
independtly on your computer. The localisation of thee tools have to be 
specified in the configuration file ``PATH.ini``. The ``PATH.ini`` file is a 
basic configuration file used by the python `ConfigParser module 
<https://docs.python.org/2/library/configparser.html>`_ .

    | ``PATH.ini`` syntax:
    | [annotation]
    | pfamscan:/opt/global/bin/pfam_scan.pl
    | [similarity]
    | ...

The following executable should be installed on your computer and present in 
the ``PATH.ini`` file:
    
    * pfam_scan.pl, script from the `Pfam ftp <ftp://ftp.sanger.ac.uk/pub/databases/Pfam/Tools/>`_
    * proteinortho, script, from the `proteinortho web site <https://www.bioinf.uni-leipzig.de/Software/proteinortho/proteinortho_v4.26.tar.gz>`_
    * hmmpress, binary from the `HMMER package <http://hmmer.janelia.org/software>`_

The Pfam database version 27.0 should also be `present <ftp://ftp.sanger.ac.uk/pub/databases/Pfam/>`_.


Install
-------

porthoDom is composed of two parts, a C program computing the pairwise domain
similarity and a python gluing everything together.

For a local installation

.. code-block:: sh

    python install setup.py --user

Be sure that the local installation path is in our python environment path
add  ``$HOME/.local/lib`` in ``PYTHONPATH``


.. code-block:: sh

    export PYTHONPATH=$HOME/.local/lib/

add the previous line to your bashrc or personal bashrc

.. code-block:: sh

    echo "export PYTHONPATH=$HOME/.local/lib/" >> ~/.bashrc

move the local compute_similarity program from bin to a directory where it will
be accessible by the system

If administrator right are provided

.. code-block:: sh

    sudo python install setup.py


