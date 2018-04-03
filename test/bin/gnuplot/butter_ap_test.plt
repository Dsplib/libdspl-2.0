
set logscale x
unset key
set grid
set xlabel "frequency, rad/s"
set ylabel "Butterworth filter magnitude, dB"
set yrange [-100:5]
plot 'dat/butter_ap_test_mag.txt' with lines
