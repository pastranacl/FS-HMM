#!/bin/bash
# Cesar Lopez Pastrana 2016

clear

#bash 
column1=\$5
column2=\$6
column3=\$2

gnuplot<< TOEND

file_name='../states.dat'

#======================================================
# COLORS
#======================================================

darkgray="#696969"
darkred="#8B0000"
lightgray="#A9A9A9"

#======================================================
# LINE STYLE
#======================================================

# For the plots
set style line 1 linetype 1 linewidth 2 linecolor rgb lightgray
set style line 2 linetype 1 linewidth 2 linecolor rgb darkred 
set style line 3 linetype 1 linewidth 2 linecolor rgb "blue"
#set style line 4 linecolor rgb lightgray pt 7   # circle

# For the axes borders
set style line 99 linetype 1 linewidth 1.5 linecolor rgb "black"

#======================================================
# FORMAT AXES
#======================================================

#set size 0.85, 1

set xlabel "Time (ms)" font " ,20
set ylabel "Extension (nm)" font " ,20"

# Axes
set yrange [1000:2000]
set xrange [0:2000]
set ytics 200 font " ,14"
set xtics 400 font " ,14" # Step size in the xticks and font size


# Apply the thickness to the border
set border linestyle 99

#======================================================
# PLOT DATA
#======================================================

###### eps ######
#set terminal postscript eps size 2,1  enhanced color font 'Arial,20' linewidth 2
#set output 'plot.eps'

###### eps LATEX ######
#set terminal epslatex size 7,3 standalone color colortext 10
#set output 'terminal_epslatex.tex'

# eps LAtex without latex document
set terminal epslatex size 7.0in,5in  font 'phv,14' standalone color lw 3 # phv sets Helvetica font
set output 'Figure1.tex'
# Then run in bash
#latex Figure1.tex
#pdflatex Figure1.tex

####### pdf ######
#set terminal postscript enhanced color size 10,5  font 'Times,28' lw 4
#set output '| ps2pdf - plot.pdf'
#set encoding iso_8859_1

#
set key right bottom

plot file_name using 1:($column1) title "" with lines ls 1,\
     file_name using 1:($column2) title "HMM" with lines ls 2,\
     file_name using 1:($column3) title "Steps" with lines ls 3
   
TOEND

# Return to bash

if  [ -e Figure1.tex ]
then
  latex Figure1.tex
  pdflatex Figure1.tex Figure1.pdf
fi

for figure_files in  `ls | grep Figure1 | grep -v pdf`
do
  rm $figure_files
done
