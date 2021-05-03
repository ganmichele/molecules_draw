#!/usr/bin/env bash

./rotate.py example.xyz -R 0 30 30 --order xzy -o example.xyz --tikz test.tex --depth --angles H1C0H2 H2C0H3 --asize 0.5 --connect connectivity_dist.txt
pdflatex test.tex 1> /dev/null

