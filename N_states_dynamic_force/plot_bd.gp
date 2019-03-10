#!/usr/bin/gnuplot -persist
# Cesar Lopez Pastrana 2016


name='bd.dat'


##### LINES ###

set style line 1 linetype 1 linecolor rgb "#555555"
set style line 2 linetype 1 linewidth 2 linecolor rgb "black"
set style line 3 linetype 1 linewidth 2 linecolor rgb "red"

# For the axes borders
set style line 100 linetype 1 linewidth 1.5 linecolor rgb "black"


####################################

# Creates the data for the histogram

set output 'hist.dat'
n = 500 #number of intervals
max = 5000  #max value
min = 0     #min value
width = (max-min)/n #interval width
set boxwidth width

#function used to map a value to the intervals
hist(x, width) = width*floor(x/width)+width/2.0

set table      
plot name using (hist($5,width)):(1.0) smooth freq w p    

unset table

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


plot name using 1:6 title "" with lines ls 1,\
     name using 1:2 title "" with lines ls 3


# Plot-2 features (the histogram)
#############################
set size 0.15, 1
set origin 0.85, 0 # position
set border linestyle 100
set xlabel "    "
set ylabel ""

stats 'hist.dat' using 1
#show variables all
set yrange [0:2000]
set xrange[10:3000]


set xtics 1000 # steps in the xticks
set format y "" 
set format x " " 
#set xtic rotate by 0
unset ytics
#set y2tics rotate by 0


plot 'hist.dat' using 2:1  title "" with line ls 2 #smooth bezier
     
  
unset multiplot

#
