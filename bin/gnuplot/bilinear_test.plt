
if(!exists("plotterm")) plotterm = 'wxt'
if(plotterm eq "pngcairo") set output 'img/bilinear.png'

unset key
set grid
set xlabel "normalized frequency"

set terminal plotterm size 520, 380 enhanced font 'Verdana,8'


set ylabel "Magnitude, dB"
set yrange [-80:5]
plot 'dat/bilinear.txt' with lines



