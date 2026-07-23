set term post eps
set output 'plot.eps'
set multiplot
set yrange [-0.025:0.025]
set xrange [20:60]
unset key
set xlabel "time"
set ylabel "amplitude"
plot "totalPointCarre.dat" using 1:2 w l lt 1
plot "totalPointCarre.dat" using 1:3 w l lt 2
plot "totalPointCarre.dat" using 1:4 w l lt 3
unset multiplot
set term x11
set output
