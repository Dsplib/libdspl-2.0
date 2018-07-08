set logscale x
unset key
set grid
set xlabel "frequency, rad/s"

set lmargin at screen 0.10

set terminal pngcairo size 560,420 enhanced font 'Verdana,8'
set output 'img/cheby1_ap_test_mag.png'
set ylabel "Chebyshev type 1 filter magnitude, dB"
set yrange [-100:5]
plot 'dat/cheby1_ap_test_mag.txt' with lines


set terminal pngcairo size 560,420 enhanced font 'Verdana,8'
set output 'img/cheby1_ap_test_phi.png'

set ylabel "Chebyshev type 1 filter phase response, rad"
unset yrange
plot 'dat/cheby1_ap_test_phi.txt' with lines


set terminal pngcairo size 560,420 enhanced font 'Verdana,8'
set output 'img/cheby1_ap_test_tau.png'

set ylabel "Chebyshev type 1 filter groupdelay, sec"
unset yrange
plot 'dat/cheby1_ap_test_tau.txt' with lines
