#!/usr/bin/gnuplot -persist


set xrange [0:11.5]
plot "test_mass.dat" using 1:2 with lines title 'mass',\
