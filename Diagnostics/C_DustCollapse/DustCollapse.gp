set term post eps color
set output 'DustCollapse.eps'

set key top right

set pointsize 2;

set xlabel "time (s)";
set ylabel "radius (cm)"

plot 'analytic.txt'   using 1:2 w l title "Exact Solution",\
     'results_1d.txt' using 1:2 pt 5 title "1D Solution", \
     'results_2d.txt' using 1:2 pt 9 title "2D Solution", \
     'results_3d.txt' using 1:2 pt 7 ps 1.0 title "3D Solution"
