set grid
set xlabel "x"

set lmargin at screen 0.10
set key left top
set terminal pngcairo size 560,380 enhanced font 'Verdana,8'
set output 'img/cheby_poly2.png'
set ylabel "C_N(x)"
set yrange [-3.5:3.5]
plot 'dat/cheby_poly2_ord1.txt' with lines, \
     'dat/cheby_poly2_ord2.txt' with lines, \
     'dat/cheby_poly2_ord3.txt' with lines, \
     'dat/cheby_poly2_ord4.txt' with lines


