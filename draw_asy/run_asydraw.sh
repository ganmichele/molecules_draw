#!/bin/bash 
#
./asydraw.py paracetamol_eq.xyz -T connect.csv --sizelab --whichatomlabel 0 1 2 5 -O test.tex --scalefact 1 --sph_s 0.45 --cyl_s 0.05 
./asytex.sh test.tex
