# plot the Sod results.  Here we assume that we have a file called sedov.out
# with columns containing x, rho, xmom, pres

set term pngcairo size 800, 1000
set output 'hse.png'

set multiplot;

set size 1, 0.33;

#set xrange [0:0.5]
set xlabel "x";

set origin 0.0, 0.666;
set ylabel "density";
plot 'bubble_256_plt02667.slice' using 1:2 notitle w l

set origin 0.0, 0.333;
set ylabel "velocity";
plot 'bubble_256_plt02667.slice' using 1:4 notitle w l

set origin 0.0, 0.0;
set ylabel "temperature";
plot 'bubble_256_plt02667.slice' using 1:3 notitle w l

unset multiplot;
set term x11;
