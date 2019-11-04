if(!exists("plotterm")) plotterm = 'wxt'
if(plotterm eq "pngcairo") set output 'img/sine_int.png'


set grid
set xlabel "x"

set lmargin at screen 0.10
set key left top
set terminal plotterm size 560,280 enhanced font 'Verdana,8'

set ylabel "Si(x), sinc(x)"
set yrange [-2:2]
plot 'dat/dat0.txt' with lines  title "Si(x)", \
'dat/dat1.txt' with lines title "sinc(x)"

