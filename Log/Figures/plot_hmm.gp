#!/usr/bin/gnuplot -persist
# Cesar Lopez Pastrana 2016

file_name = 'states.dat'

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

# For the axes borders
set style line 3 linetype 1 linewidth 1.5 linecolor rgb "black"

#======================================================
# FORMAT AXES
#======================================================

set size 0.85, 1

set xlabel "Time (ms)" font " ,18
set ylabel "Extension ( {/Symbol m}m)" font " ,18"

# Axes
set yrange [1:2]

set ytics 0.3 font " ,14"
set xtics 1000 font " ,14" # Step size in the xticks and font size


# Apply the thickness to the border
set border linestyle 3

#======================================================
# PLOT DATA
#======================================================

###### eps ######
#set terminal postscript eps size 2,1  enhanced color font 'Arial,20' linewidth 2
#set output 'plot.eps'

###### eps LATEX ######
set terminal epslatex size 3.5,2.62 standalone color colortext 10
set output 'terminal_epslatex.tex'

# eps LAtex without latex document
set terminal epslatex size 11inch,9inch standalone color colortext 10
set output 'Figure1.tex'
# Then run in bash
#latex Figure1.tex
#pdflatex Figure1.tex

####### pdf ######
#set terminal postscript enhanced color size 10,5  font 'Times,28' lw 4
#set output '| ps2pdf - plot.pdf'
#set encoding iso_8859_1

plot file_name using 1:($5/1000) title "" with lines ls 1,\
     file_name using 1:($6/1000) title "" with lines ls 2