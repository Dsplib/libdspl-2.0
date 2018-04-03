
set logscale x
unset key
set grid
set xlabel "frequency, rad/s"
set ylabel "Chebyshev type 1 filter magnitude, dB"
set yrange [-100:5]
plot 'dat/cheby1_ap_test_mag.txt' with lines
