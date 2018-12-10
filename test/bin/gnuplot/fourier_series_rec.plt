### Start multiplot (2x2 layout)

set terminal wxt 0 size 800,640 enhanced font 'Verdana,8' position 0,0
unset key
set multiplot layout 4,2 rowsfirst
#
plot 'dat/fourier_series_pimp0.txt' with lines,\
	'dat/fourier_series_pimp_rec_5.txt' with lines
#
plot 'dat/fourier_series_saw0.txt' with lines,\
	'dat/fourier_series_saw_rec_5.txt' with lines
#
plot 'dat/fourier_series_pimp0.txt' with lines,\
	'dat/fourier_series_pimp_rec_9.txt' with lines

#
plot 'dat/fourier_series_saw0.txt' with lines,\
	'dat/fourier_series_saw_rec_9.txt' with lines
#
plot 'dat/fourier_series_pimp0.txt' with lines,\
	'dat/fourier_series_pimp_rec_21.txt' with lines

#
plot 'dat/fourier_series_saw0.txt' with lines,\
	'dat/fourier_series_saw_rec_21.txt' with lines
#
plot 'dat/fourier_series_pimp0.txt' with lines,\
	'dat/fourier_series_pimp_rec_61.txt' with lines

#
plot 'dat/fourier_series_saw0.txt' with lines,\
	'dat/fourier_series_saw_rec_61.txt' with lines

unset multiplot




set terminal pngcairo size 800,640 enhanced font 'Verdana,8'
set output 'img/fourier_series_rec.png'
unset key
set multiplot layout 4,2 rowsfirst
#
plot 'dat/fourier_series_pimp0.txt' with lines,\
	'dat/fourier_series_pimp_rec_5.txt' with lines
#
plot 'dat/fourier_series_saw0.txt' with lines,\
	'dat/fourier_series_saw_rec_5.txt' with lines
#
plot 'dat/fourier_series_pimp0.txt' with lines,\
	'dat/fourier_series_pimp_rec_9.txt' with lines

#
plot 'dat/fourier_series_saw0.txt' with lines,\
	'dat/fourier_series_saw_rec_9.txt' with lines
#
plot 'dat/fourier_series_pimp0.txt' with lines,\
	'dat/fourier_series_pimp_rec_21.txt' with lines

#
plot 'dat/fourier_series_saw0.txt' with lines,\
	'dat/fourier_series_saw_rec_21.txt' with lines
#
plot 'dat/fourier_series_pimp0.txt' with lines,\
	'dat/fourier_series_pimp_rec_61.txt' with lines

#
plot 'dat/fourier_series_saw0.txt' with lines,\
	'dat/fourier_series_saw_rec_61.txt' with lines
unset multiplot
