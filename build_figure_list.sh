#!/bin/bash


grep \includegraphics ./latex/*tex | grep -v "%" | cut -f 2 -d"{" | cut -f1 -d"}" > figure_list.txt
