
if(!exists("plotterm")) plotterm = 'wxt'
if(plotterm eq "pngcairo") set output 'img/gnuplot_script.png'

set grid

set xlabel "x"

set terminal plotterm size 560,280 enhanced font 'Verdana,8'

set ylabel "sin(x), cos(x)"

plot 'dat/sine.txt' with lines title "sin(x)", \
     'dat/cosine.txt' with lines  title "cos(x)"

