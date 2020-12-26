#!/usr/bin/gnuplot -persist

set xrange [0:3]
set yrange [0:3]
set xlabel "абсцисса"
set ylabel "ордината"
# set pm3d map
 
# set palette rgbformulae 30,31,32
# set cbrange [0.99:1.01]
plot "test2.dat" every 2:2:0:0 using 1:2:($3/4):($4/4)\
with vectors title 'velocity'

