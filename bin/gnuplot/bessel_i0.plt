set grid
set xlabel "x"

set lmargin at screen 0.10
set key left top
set terminal pngcairo size 560,380 enhanced font 'Verdana,8'
set output 'img/bessel_i0.png'
set ylabel "I_0(x)"
set yrange [0:5]
plot 'dat/dat0.txt' with lines

