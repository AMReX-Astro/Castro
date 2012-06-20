# plot the Sod results.  Here we assume that we have files 
# called sodx.out, sody.out, and sodz.out

set term post eps color
set output 'sod_3d.eps'

set multiplot;

set size 0.5, 0.5;

set xlabel "x/y/z";

set style line 1  lw 5

set origin 0.0, 0.5;
set key top right
set ylabel "density";
plot 'sodx.out' using 1:2 title 'x' with points 1 1,\
     'sody.out' using 1:2 title 'y' with points 2 2,\
     'sodz.out' using 1:2 title 'z' with points 3 4,\
     'sod-exact.out' using 1:2 notitle with lines lt 1 lw 2;

set origin 0.5, 0.5;
set key top left
set ylabel "velocity";
plot 'sodx.out' using 1:($3/$2) title 'x' with points 1 1,\
     'sody.out' using 1:($4/$2) title 'y' with points 2 2,\
     'sodz.out' using 1:($5/$2) title 'z' with points 3 4,\
     'sod-exact.out' using 1:3 notitle with lines lt 1 lw 2;


set origin 0.0, 0.0;
set key top right
set ylabel "pressure";
plot 'sodx.out' using 1:10 title 'x' with points 1 1,\
     'sody.out' using 1:10 title 'y' with points 2 2,\
     'sodz.out' using 1:10 title 'z' with points 3 4,\
     'sod-exact.out' using 1:4 notitle with lines lt 1 lw 2;

set origin 0.5, 0.0;
set key top left
set ylabel "internal energy";
plot 'sodx.out' using 1:($7/$2) title 'x' with points 1 1,\
     'sody.out' using 1:($7/$2) title 'y' with points 2 2,\
     'sodz.out' using 1:($7/$2) title 'z' with points 3 4,\
     'sod-exact.out' using 1:5 notitle with lines lt 1 lw 2;

unset multiplot;
set term x11;
