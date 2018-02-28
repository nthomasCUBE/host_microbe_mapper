#!/bin/bash

PATH=/naslx/projects_mpiio/pr74ma/ge34juq2/apps/reago-1.1-release-2015.12.18/gt-1.5.10-Linux_x86_64-64bit/bin/:$PATH
PATH=/naslx/projects_mpiio/pr74ma/ge34juq2/apps/reago-1.1-release-2015.12.18/infernal-1.1.2-linux-intel-gcc/binaries/:$PATH

python filter_input.py sample_1.fasta sample_2.fasta filter_out cm ba 10
python reago.py filter_out/filtered.fasta xxx -l 101



