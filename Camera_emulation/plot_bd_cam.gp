#!/usr/bin/gnuplot -persist
# Cesar Lopez Pastrana 2016

bd_raw='bd.dat'
bd_cam='bd_cam.dat'

##### LINES ###

set style line 1 linetype 1 linecolor rgb "#555555"
set style line 2 linetype 1 linewidth 2 linecolor rgb "black"
set style line 3 linetype 1 linewidth 2 linecolor rgb "red"

# For the axes borders
set style line 100 linetype 1 linewidth 1.5 linecolor rgb "black"


#####################################

# Comment these lines to avoid saving and thus show
# the plot after execution

# For SVG
#set terminal svg size 350,262 fname 'FreeSerif' fsize 12
#set output 'bd.svg'

# For PNG
set terminal png transparent font FreeSerif 18 size 800,600
set output 'bd.png'

set multiplot layout 1, 2

###################################
# Plot-1 features (the dynamics)
###################################
set size 0.88, 1
set border linestyle 100
set xlabel "Time (ms)" font " ,14"
set ylabel "Extension (nm)" font " ,14"

set yrange [0:2000]
set ytics 300 font " ,12"
set xtics 1000 font " ,12" # steps in the xticks


plot bd_raw using 1:5 title "" with lines ls 1,\
     bd_cam using 1:5 title "" with lines ls 99
    

#
