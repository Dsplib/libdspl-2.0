
set terminal wxt 0 size 560,480 enhanced font 'Verdana,8' position 0,0
unset key
set grid
set lmargin 8
set multiplot layout 2,1 rowsfirst
set xlabel 'Частота, рад/с'
#
set ylabel 'Амплитудный спектр'
plot[-10*pi:10*pi] 'dat/fourier_series_pimp_mag.txt' with  impulses lt 1 ,\
		'dat/fourier_series_pimp_mag.txt' with points pt 7 ps 0.5 lt 1
#
set ylabel 'Фазовый спектр'
plot[-10*pi:10*pi] 'dat/fourier_series_pimp_phi.txt' with  impulses lt 1 ,\
	'dat/fourier_series_pimp_phi.txt' with points pt 7 ps 0.5 lt 1


unset multiplot


set terminal pngcairo size 560,480 enhanced font 'Verdana,8'
set output 'img/fourier_series_pimp.png'
set multiplot layout 2,1 rowsfirst
#
plot[-10*pi:10*pi] 'dat/fourier_series_pimp_mag.txt' with  impulses lt 1 ,\
	'dat/fourier_series_pimp_mag.txt' with points pt 7 ps 0.5 lt 1
#
plot[-10*pi:10*pi] 'dat/fourier_series_pimp_phi.txt' with  impulses lt 1 ,\
	'dat/fourier_series_pimp_phi.txt' with points pt 7 ps 0.5 lt 1

unset multiplot
