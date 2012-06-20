# plot the Sod results.  Here we assume that we have a file called sedov.out
# with columns containing x, rho, xmom, pres

set term post eps color
set output 'sedov.eps'

set multiplot;

set size 1, 0.33;

set xrange [0:0.5]
set xlabel "r";

set origin 0.0, 0.666;
set ylabel "density";
plot 'sedov.out' using 1:2 notitle,\
     'spherical_sedov.dat' using 2:3 notitle with lines;

set origin 0.0, 0.333;
set ylabel "velocity";
plot 'sedov.out' using 1:3 notitle,\
     'spherical_sedov.dat' using 2:6 notitle with lines;

set origin 0.0, 0.0;
set ylabel "pressure";
plot 'sedov.out' using 1:4 notitle,\
     'spherical_sedov.dat' using 2:5 notitle with lines;

unset multiplot;
set term x11;
