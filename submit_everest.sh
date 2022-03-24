#!/bin/bash
#$ -S /bin/bash
#$ -q 12cores.q,36cores.q
#$ -cwd
#$ -o tmp/output
#$ -e tmp/error
#$ -j n
#$ -r y
#$ -N InterfaceExtension
conda deactivate
conda activate
python main.py -pdb 6waq_fixed.pdb -qc A -ic B -sr 4.4 -mi 1 -r result.txt


