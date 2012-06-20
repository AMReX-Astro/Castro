# plot the Test3 results.  Here we assume that we have files 
# called test3x.out, test3y.out, and test3z.out

set term post eps color
set output 'test3_3d.eps'

set multiplot;

set size 0.5, 0.5;

set xlabel "x/y/z";

set style line 1  lw 5

set origin 0.0, 0.5;
set key top left
set ylabel "density";
plot 'test3x.out' using 1:2 title 'x' with points 1 1,\
     'test3y.out' using 1:2 title 'y' with points 2 2,\
     'test3z.out' using 1:2 title 'z' with points 3 4,\
     'test3-exact.out' using 1:2 notitle with lines lt 1 lw 2;

set origin 0.5, 0.5;
set key top left
set ylabel "velocity";
plot 'test3x.out' using 1:($3/$2) title 'x' with points 1 1,\
     'test3y.out' using 1:($4/$2) title 'y' with points 2 2,\
     'test3z.out' using 1:($5/$2) title 'z' with points 3 4,\
     'test3-exact.out' using 1:3 notitle with lines lt 1 lw 2;

set origin 0.0, 0.0;
set key top right
set ylabel "pressure";
plot 'test3x.out' using 1:10 title 'x' with points 1 1,\
     'test3y.out' using 1:10 title 'y' with points 2 2,\
     'test3z.out' using 1:10 title 'z' with points 3 4,\
     'test3-exact.out' using 1:4 notitle with lines lt 1 lw 2;

set origin 0.5, 0.0;
set key top right
set ylabel "internal energy";
plot 'test3x.out' using 1:($7/$2) title 'x' with points 1 1,\
     'test3y.out' using 1:($7/$2) title 'y' with points 2 2,\
     'test3z.out' using 1:($7/$2) title 'z' with points 3 4,\
     'test3-exact.out' using 1:5 notitle with lines lt 1 lw 2;

unset multiplot;
set term x11;
