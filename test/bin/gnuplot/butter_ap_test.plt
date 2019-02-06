set logscale x
unset key
set grid
set xlabel "frequency, rad/s"

set terminal pngcairo size 920, 260 enhanced font 'Verdana,8'
set output 'img/butter_ap_test.png'


set multiplot layout 1,3 rowsfirst
set ylabel "Magnitude, dB"
set yrange [-100:5]
plot 'dat/butter_ap_test_mag.txt' with lines

set ylabel "Phase response, rad"
unset yrange
plot 'dat/butter_ap_test_phi.txt' with lines

set ylabel "Groupdelay, sec"
unset yrange
plot 'dat/butter_ap_test_tau.txt' with lines

unset multiplot

