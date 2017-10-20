#!/usr/bin/gnuplot -persist

set grid

#set logscale xy

set xlabel "VMC_{variance}"
set ylabel "DMC_{variance}"

f(x) = a*x
fit f(x) "vmc_dmc_variance.dat" using 2:4 via a

stats "vmc_dmc_variance.dat" using 2:4 name "fit"

plot "vmc_dmc_variance.dat" using 2:4:3:5 with xyerrorbars notitle, f(x) notitle, \
     f(x) title sprintf("y = %.3f*x (r = %.3f)", a, fit_correlation), \
     "vmc_dmc_variance.dat" using 2:4:1 with labels offset -2.0,-0.5 font "arial,8" notitle
