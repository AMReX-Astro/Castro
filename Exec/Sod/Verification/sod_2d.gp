# plot the Sod results.  Here we assume that we have files 
# called sodx.out and sody.out

set term post eps color
set output 'sod_2d.eps'

set multiplot;

set size 1, 0.33;

set xlabel "x/y";

set style line 1  lw 5

set origin 0.0, 0.666;
set ylabel "density";
plot 'sodx.out' using 1:2 title 'x' with points 1 1,\
     'sody.out' using 1:2 title 'y' with points 2 2,\
     'sod-exact.out' using 1:2 notitle with lines ls 1;

set origin 0.0, 0.333;
set ylabel "velocity";
plot 'sodx.out' using 1:($3/$2) title 'x' with points 1 1,\
     'sody.out' using 1:($4/$2) title 'y' with points 2 2,\
     'sod-exact.out' using 1:3 notitle with lines ls 1;

set origin 0.0, 0.0;
set ylabel "pressure";
plot 'sodx.out' using 1:9 title 'x' with points 1 1,\
     'sody.out' using 1:9 title 'y' with points 2 2,\
     'sod-exact.out' using 1:4 notitle with lines ls 1;

unset multiplot;
set term x11;
