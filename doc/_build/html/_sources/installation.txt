.. _installation


Requirement
===========


the following softwares and database must be installed

pfam_scan.pl
hmm database pfam-27.0
hmmscan 
proteinortho


Install
=======


porthoDom is composed of two parts, a C program computing the pairwise domain
similarity and a python gluing everything together.

For a local installation

python install setup.py --user

Be sure that the local installation path is in our python environment path
add  $HOME/.local/lib in PYTHONPATH

export PYTHONPATH=$HOME/.local/lib/

add the previous line to your bashrc or personal bashrc

echo "export PYTHONPATH=$HOME/.local/lib/" >> ~/.bashrc

move the local compute_similarity program from bin to a directory where it will
be accessible by the system

If administrator right are provided

sudo python install setup.py


