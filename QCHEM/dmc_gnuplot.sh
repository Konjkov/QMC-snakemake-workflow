#!/usr/bin/gnuplot -persist

f(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12) = h*x1 + be*x2 + b*x3 + c*x4 + n*x5 + o*x6 + f*x7 + al*x8 + si*x9 + p*x10 + s*x11 + cl*x12

# Initial values
h  =  -0.5
be = -14.66736
b  = -24.65393
c  = -37.8450
n  = -54.5893
o  = -75.067
f  = -99.734
al = -242.3460
si = -289.3590
p  = -341.2590
s  = -398.1100
cl = -460.1480

# Fitting
set dummy x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12
fit f(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12) "dmc_energy.dat" using 2:3:4:5:6:7:8:9:10:11:12:13:14 via h, be, b, c, n, o, f, al, si, p, s, cl

Min = -10.0 # where binning starts
Max =  35.0 # where binning ends
n = 18 # the number of bins
binwidth = (Max-Min)/n # binwidth
bin(x) = binwidth*(floor((x-Min)/binwidth)+0.5) + Min

set xlabel "TAE(literature)-TAE(calculated)(kcal/mol)"
set ylabel "molecules"

set boxwidth binwidth

#set term png
#set output "dmc_energy.png"

plot 'dmc_energy.dat' using (bin($16)):(1.0) smooth freq with boxes notitle
