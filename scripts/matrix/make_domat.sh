#!/bin/bash
# collection of script to create the domain similarity matrix
# WARNING, hhsearch only works with hhr formated hmms
# Pfam-A.hhr models can be downloaded at
# http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/
# for other database the hhr models need to be generated using hhmake

hmmdb=$1                 # path to the file with all the hmmdata, example Pfam-A.hmm
workout=$2               # path to the working directory
matrix=domsim_matrix.dat # output name
simcutoff=1              # simlarity cutoff

# bunch of directory and subdirectories to create
if [ ! -d "$workout" ]; then
  mkdir -p $workout
fi

tmphmm=$workout/tmphmm
tmphhsearch=$workout/tmphhsearch

if [ ! -d "$tmphmm" ]; then
  mkdir -p $tmphmm
fi

if [ ! -d "$tmphhsearch" ]; then
  mkdir -p $tmphhsearch
fi

## important stuff starts here ##

# cut the inital pfam file into smaller pieces
python cut_seed.py -i $hmmdb -o $tmphmm -m ${workout}/map_name2id.dat

# run for each pfam pieces hhsearch against the hmmdatabase
for hmm in `ls $tmphmm`
do
    hhsearch -i $tmphmm/$hmm -d $hmmdb -scores $tmphhsearch/${hmm%.*}.scores &> $tmphhsearch/${hmm%.*}.log
done

# read hhsearch scores and make domain matrix
python make_matrix.py -i $tmphhsearch/*.scores -o $matrix -n matbin_$simcutoff -s PROBAB -m ${workout}/map_name2id.dat -t $simcutoff

# WARNING if this script is not working
# python make_matrix2.py -i $tmphhsearch/*.scores -o $workout/tmp_matrix.dat -d PROBAB -m ${workout}/map_name2id.dat 
# mat2bincrs $workout/tmp_matrix.dat $matrix $simcutoff matbin_$simcutoff

    
