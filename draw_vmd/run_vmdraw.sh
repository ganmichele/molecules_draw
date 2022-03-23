#!/bin/bash

if [ -z $1 ]; then
    input="vmd_commands.vmd"
else
    input="$1"
fi

./vmdraw.py paracetamol_eq.xyz -T connect_fat.csv --displ no -O $input  

vmd -e "$input"
