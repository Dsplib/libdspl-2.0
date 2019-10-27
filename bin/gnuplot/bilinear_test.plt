unset key
set grid
set xlabel "normalized frequency"

set terminal pngcairo size 520, 380 enhanced font 'Verdana,8'
set output 'img/bilinear.png'

set ylabel "Magnitude, dB"
set yrange [-80:5]
plot 'dat/bilinear.txt' with lines



