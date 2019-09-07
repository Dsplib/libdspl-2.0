unset key
set grid
set xlabel " normalized frequency"

set terminal pngcairo size 920, 840 enhanced font 'Verdana,8'
set output 'img/iir_test.png'
set ylabel "Magnitude, dB"
set yrange [-100:5]

set multiplot layout 4,4 rowsfirst

plot 'dat/iir_butter_lpf.txt' with lines
plot 'dat/iir_butter_hpf.txt' with lines
plot 'dat/iir_butter_bpf.txt' with lines
plot 'dat/iir_butter_bsf.txt' with lines

plot 'dat/iir_cheby1_lpf.txt' with lines
plot 'dat/iir_cheby1_hpf.txt' with lines
plot 'dat/iir_cheby1_bpf.txt' with lines
plot 'dat/iir_cheby1_bsf.txt' with lines

plot 'dat/iir_cheby2_lpf.txt' with lines
plot 'dat/iir_cheby2_hpf.txt' with lines
plot 'dat/iir_cheby2_bpf.txt' with lines
plot 'dat/iir_cheby2_bsf.txt' with lines

plot 'dat/iir_ellip_lpf.txt' with lines
plot 'dat/iir_ellip_hpf.txt' with lines
plot 'dat/iir_ellip_bpf.txt' with lines
plot 'dat/iir_ellip_bsf.txt' with lines

unset multiplot