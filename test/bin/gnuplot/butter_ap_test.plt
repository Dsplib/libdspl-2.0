set logscale x
unset key
set grid
set xlabel "frequency, rad/s"

set lmargin at screen 0.10


set terminal wxt 0 size 460,320 enhanced font 'Verdana,8'

set ylabel "Butterworth filter magnitude, dB"
set yrange [-100:5]
plot 'dat/butter_ap_test_mag.txt' with lines


set terminal pngcairo size 560,420 enhanced font 'Verdana,8'
set output 'img/butter_ap_test_mag.png'
replot


set terminal wxt 0 size 460,320 enhanced font 'Verdana,8'
set ylabel "Butterworth filter phase response, rad"
unset yrange
plot 'dat/butter_ap_test_phi.txt' with lines

set terminal pngcairo size 560,420 enhanced font 'Verdana,8'
set output 'img/butter_ap_test_phi.png'
replot



set terminal wxt 0 size 460,320 enhanced font 'Verdana,8'
set ylabel "Butterworth filter groupdelay, sec"
unset yrange
plot 'dat/butter_ap_test_tau.txt' with lines

set terminal pngcairo size 560,420 enhanced font 'Verdana,8'
set output 'img/butter_ap_test_tau.png'
replot
