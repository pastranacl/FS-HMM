#!/bin/bash

clear 
LANG=C

file=bd.dat
lines=`cat $file | wc -l`
if [ -e ./cam_em ]
then
  rm ./cam_em
fi
gcc -std=c99 ./cam_em.c -o cam_em
./cam_em $lines
#./cam_em.c $lines
./plot_bd_cam.gp