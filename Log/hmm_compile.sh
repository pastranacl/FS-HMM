clear

if [ -e ./states.dat ]; then
  rm ./viterbi.dat
  rm ./states.dat
fi

gcc -O2 -lm -std=c99 hmm.c energy.c files.c -o hmm
START=$(date +%s)
./hmm
END=$(date +%s)
DIFF=$(( $END - $START ))
echo -e "\n"
echo "Execution time: $DIFF s"


paste -d' ' bd.dat viterbi.dat >> states.dat


./plot_hmm.gp

# To do the DLL
#gcc -c -fpic -O3 -lm -std=c99 hmm.c
#gcc -shared -o libhmm.so hmm.o