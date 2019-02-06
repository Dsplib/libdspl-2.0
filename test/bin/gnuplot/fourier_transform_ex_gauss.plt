set terminal pngcairo size 560,480 enhanced font 'Verdana,8'
set output 'img/fourier_series_pimp.png'
unset key
set grid
set lmargin 8
set multiplot layout 2,1 rowsfirst

#
set xlabel 't, c'
set ylabel 's(t)'
plot[-10:10] 'dat/fourier_transform_ex_gauss_time_0.5.txt' with lines,\
		     'dat/fourier_transform_ex_gauss_time_1.0.txt' with lines,\
             'dat/fourier_transform_ex_gauss_time_2.0.txt' with lines
#
set xlabel 'w, рад/c'
set ylabel 'S(w)'
plot[-4*pi:4*pi] 'dat/fourier_transform_ex_gauss_freq_0.5.txt' with lines,\
		         'dat/fourier_transform_ex_gauss_freq_1.0.txt' with lines,\
                 'dat/fourier_transform_ex_gauss_freq_2.0.txt' with lines


unset multiplot


