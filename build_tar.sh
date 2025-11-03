#!/bin/bash

NAME="varEstRisk2025lawless"

echo "making copy of core files"
mkdir temp_build
# cp latex/head.tex temp_build/
# cp -r latex/images temp_build/
cp latex/references.bib temp_build/
cp latex/${NAME}.pdf temp_build/
cp latex/${NAME}.tex temp_build/

echo "zipping tar.gz"
tar -czvf ${NAME}.tar.gz temp_build

echo "removing temp_build"
rm -r temp_build
