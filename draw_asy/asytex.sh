#!/bin/bash

doc="program to draw pictures with latex & asy"
usage="$0 file.tex"

if [ -z $1 ]; then
    echo $doc
    echo "you must pass exactly one argument file"
    echo $usage
fi

fname=`echo $1 | cut -d'.' -f1`
texname="${fname}.tex"
asyname="${fname}-1.asy"

#echo "$fname"
#echo "$texname"
#echo "$asyname"

pdflatex "$texname" 1> /dev/null &&
echo "pdflatex 1 done"
asy "$asyname" &&
echo "asy done"
pdflatex "$texname" 1> /dev/null
echo "pdflatex 2 done"
