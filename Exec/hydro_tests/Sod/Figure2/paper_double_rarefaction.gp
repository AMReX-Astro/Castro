# plot the Test2 results.  Here we assume that we have files 
# called test2x.out, test2y.out, and test2z.out

set term post eps color
set output 'figure2.eps'

set multiplot;

set size 0.5, 0.5;

set xlabel "x";

set style line 1 lt rgb "black" lw 1

set origin 0.0, 0.5;
set key top center
set ylabel "density";
plot '128/den_0015.pla' using 1:2 title   '128' with p pt 1 ps 1.0,\
     '512/den_0015.pla' using 1:2 title   '512' with p pt 2 ps 0.5,\
     '2048/den_0015.pla' using 1:2 title '2048' with p pt 8 ps 0.3,\
      'test2-exact.out' using 1:2 title 'exact' with lines ls 1;
set origin 0.5, 0.5;
set key top left
set ylabel "velocity";
plot '128/vel_0015.pla' using 1:2 title   '128' with p pt 1 ps 1.0,\
     '512/vel_0015.pla' using 1:2 title   '512' with p pt 2 ps 0.5,\
     '2048/vel_0015.pla' using 1:2 title '2048' with p pt 8 ps 0.3,\
      'test2-exact.out' using 1:3 title 'exact' with lines ls 1;

set origin 0.0, 0.0;
set key top center
set ylabel "pressure";
plot '128/pres_0015.pla' using 1:2 title   '128' with p pt 1 ps 1.0,\
     '512/pres_0015.pla' using 1:2 title   '512' with p pt 2 ps 0.5,\
     '2048/pres_0015.pla' using 1:2 title '2048' with p pt 8 ps 0.3,\
      'test2-exact.out' using 1:4 title 'exact' with lines ls 1;

set origin 0.5, 0.0;
set key top center
set ylabel "internal energy";
plot '128/eint_0015.pla' using 1:2 title   '128' with p pt 1 ps 1.0,\
     '512/eint_0015.pla' using 1:2 title   '512' with p pt 2 ps 0.5,\
     '2048/eint_0015.pla' using 1:2 title '2048' with p pt 8 ps 0.3,\
      'test2-exact.out' using 1:5 title 'exact' with lines ls 1;

unset multiplot;
set term x11;
