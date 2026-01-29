# plot the Sod results.  Here we assume that we have a file called sedov.out
# with columns containing x, rho, xmom, pres

set term post eps color solid enhanced
set output 'sedov_sph.eps'
set colorsequence classic

set multiplot;

set size 0.5, 0.5;
set key right top
set xrange [0:0.3]
set xlabel "r";

set origin 0.0, 0.5;
set ylabel "density";
plot 'sedov_1d_sph.out' using 1:2 title '1-d',\
     'sedov_2d_sph_in_cyl.out' using 1:2 title '2-d cyl',\
     'sedov_2d_sph_in_sph.out' using 1:2 title '2-d sph',\
     'sedov_3d_sph.out' using 1:2 title '3-d',\
     'spherical_sedov.dat' using 2:3 notitle with lines;

set origin 0.5, 0.5;
set ylabel "velocity";
plot 'sedov_1d_sph.out' using 1:3 title '1-d',\
     'sedov_2d_sph_in_cyl.out' using 1:3 title '2-d cyl',\
     'sedov_2d_sph_in_sph.out' using 1:3 title '2-d sph',\
     'sedov_3d_sph.out' using 1:3 title '3-d',\
     'spherical_sedov.dat' using 2:6 notitle with lines;

set origin 0.0, 0.0;
set ylabel "pressure";
plot 'sedov_1d_sph.out' using 1:4 title '1-d',\
     'sedov_2d_sph_in_cyl.out' using 1:4 title '2-d cyl',\
     'sedov_2d_sph_in_sph.out' using 1:4 title '2-d sph',\
     'sedov_3d_sph.out' using 1:4 title '3-d',\
     'spherical_sedov.dat' using 2:5 notitle with lines;

set origin 0.5, 0.0;
set logscale y
set format y "10^{%L}"

set ylabel "internal energy";
plot 'sedov_1d_sph.out' using 1:5 title '1-d',\
     'sedov_2d_sph_in_cyl.out' using 1:5 title '2-d cyl',\
     'sedov_2d_sph_in_sph.out' using 1:5 title '2-d sph',\
     'sedov_3d_sph.out' using 1:5 title '3-d',\
     'spherical_sedov.dat' using 2:4 notitle with lines;

unset multiplot;
set term x11;
