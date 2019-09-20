unset key
set grid
set xlabel "n"

set terminal pngcairo size 820, 340 enhanced font 'Verdana,8'
set output 'img/filter_iir_test.png'
set ylabel "s(n)"
set yrange [-3:3]

set multiplot layout 2,1 rowsfirst

plot 'dat/s.txt'  with lines

set ylabel "s_f(n)"
plot 'dat/sf.txt' with lines


unset multiplot