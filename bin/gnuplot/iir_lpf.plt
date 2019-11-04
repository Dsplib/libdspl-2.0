if(!exists("plotterm")) plotterm = 'wxt'
if(plotterm eq "pngcairo") set output 'img/iir_lpf.png'

unset key
set grid
set xlabel " normalized frequency"

set terminal plotterm size 920, 260 enhanced font 'Verdana,8'

set multiplot layout 1,3 rowsfirst
set ylabel "Magnitude, dB"
set yrange [-100:5]
plot 'dat/iir_lpf_mag.txt' with lines

set ylabel "Phase response, rad"
unset yrange
plot 'dat/iir_lpf_phi.txt' with lines

set ylabel "Groupdelay, samples"
unset yrange
plot 'dat/iir_lpf_tau.txt' with lines

unset multiplot