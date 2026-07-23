set term post eps
set output 'cross.eps'
set multiplot
set yrange [-0.000025:0.000025]
set xrange [-5:5]
unset key
set xlabel "time"
set ylabel "amplitude"
plot "CrossComputation.dat" using 1:2 w l lt 1
plot "CrossComputation.dat" using 1:4 w l lt 2
plot "CrossComputation.dat" using 1:6 w l lt 3
unset multiplot
set term x11
set output


