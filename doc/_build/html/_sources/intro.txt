.. _introduction

Introduction
============

porthoDom is a wrapper around proteinortho version 4.26.
Latter version are currently not supported but we are working on it.

porthoDom performs the domain annotation of your protein sequences or can use 
pre-existing annotation in pfam format

The domain arrangements present in the input protein sequences set are clustered
together based on our new domain similarity distance.

The clusters are used as sub-space of search for orthologous candidates, ie
proteins sequences with similar domain arrangement are grouped together per 
species and proteinortho is runned on these sub-groups.

The results of the different proteinortho runs are gathered together and a 
standard proteinortho output file is created.


