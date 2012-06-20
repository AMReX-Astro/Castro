# plot the Test2 results.  Here we assume that we a have file
# called test2x.out

set term post eps color
set output 'test2_1d.eps'

set multiplot;

set size 1, 0.33;

set xlabel "x";

set style line 1  lw 5

set origin 0.0, 0.666;
set ylabel "density";
plot 'test2x.out' using 1:2 title 'x' with points 1 1,\
     'test2-exact.out' using 1:2 notitle with lines ls 1;

set origin 0.0, 0.333;
set ylabel "velocity";
plot 'test2x.out' using 1:($3/$2) title 'x' with points 1 1,\
     'test2-exact.out' using 1:3 notitle with lines ls 1;

set origin 0.0, 0.0;
set ylabel "pressure";
plot 'test2x.out' using 1:8 title 'x' with points 1 1,\
     'test2-exact.out' using 1:4 notitle with lines ls 1;

unset multiplot;
set term x11;
