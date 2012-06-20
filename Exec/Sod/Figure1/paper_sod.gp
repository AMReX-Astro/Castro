# plot the Sod results.  

set term post eps color
set output 'figure1.eps'

set multiplot;

set size 0.5, 0.5;

set xlabel "x";

set style line 1 lt rgb "black" lw 1

set origin 0.0, 0.5;
set key top center
set ylabel "density";
plot '128/den_0034.pla' using 1:2 title   '128' with p pt 1 ps 1.0,\
     'sod-exact.out' using 1:2 title 'exact' with lines ls 1;
set origin 0.5, 0.5;
set key top left
set ylabel "velocity";
plot '128/vel_0034.pla' using 1:2 title   '128' with p pt 1 ps 1.0,\
     'sod-exact.out' using 1:3 title 'exact' with lines ls 1;

set origin 0.0, 0.0;
set key top center
set ylabel "pressure";
plot '128/pres_0034.pla' using 1:2 title   '128' with p pt 1 ps 1.0,\
     'sod-exact.out' using 1:4 title 'exact' with lines ls 1;

set origin 0.5, 0.0;
set key top center
set ylabel "internal energy";
plot '128/eint_0034.pla' using 1:2 title   '128' with p pt 1 ps 1.0,\
     'sod-exact.out' using 1:5 title 'exact' with lines ls 1;

unset multiplot;
set term x11;
