#!/usr/bin/gnuplot -persist

set grid

set xlabel "HF energy"
set ylabel "(VMC energy - HF energy)/VMC variance"

plot "hf_vmc_energy.dat" using (-$2):(($3-$2)/$4) notitle,\
     "hf_vmc_energy.dat" using (-$2):(($3-$2)/$4):1 with labels offset -2.0,-0.5 font "arial,8" notitle
