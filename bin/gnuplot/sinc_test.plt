if(!exists("plotterm")) plotterm = 'wxt'
if(plotterm eq "pngcairo") set output 'img/sinc_test.png'


set grid
set xlabel "x"

set lmargin at screen 0.10

set terminal plotterm size 560,280 enhanced font 'Verdana,8'

set ylabel "sinc(x,a)"
set yrange [-0.25:1.1]
plot 'dat/sinc_test_1.0.txt' with lines title "a = 1.0", \
'dat/sinc_test_pi.txt' with lines  title "a = pi", \
'dat/sinc_test_2pi.txt' with lines title "a = 2pi"

