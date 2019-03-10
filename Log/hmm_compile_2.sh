#!/bin/bash

# Change the languaje to english
LANG=C

clear
rm states.dat

# Compile
comp_exit=`gcc -O2 -lm -std=c99 hmm.c energy.c files.c -o hmm`

error_code=$? 

if [ error_code ]
then
	echo "Compilation completed"
	./hmm
	paste -d' ' bd.dat viterbi.dat >> states.dat
	./Figures/plot_hmm.gp
else
	echo "Compilation failed"
	echo $comp_exit
fi

# To do the DLL
#gcc -c -fpic -O3 -lm -std=c99 hmm.c
#gcc -shared -o libhmm.so hmm.o
