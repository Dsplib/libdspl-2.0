/*
* Copyright (c) 2015-2020 Sergey Bakhurin
* Digital Signal Processing Library [http://dsplib.org]
*
* This file is part of libdspl-2.0.
*
* is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser  General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* DSPL is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/



#ifdef WIN_OS
#include <windows.h>
#endif  /* WIN_OS */

#ifdef LINUX_OS
#include <dlfcn.h>
#endif /* LINUX_OS */


#include <stdio.h>
#include "dspl.h"


#ifndef BUILD_LIB

p_acos_cmplx                            acos_cmplx                    ;
p_addlog                                addlog                        ;
p_array_scale_lin                       array_scale_lin               ;
p_asin_cmplx                            asin_cmplx                    ;

p_butter_ap                             butter_ap                     ;
p_bessel_i0                             bessel_i0                     ;
p_bilinear                              bilinear                      ;
p_butter_ap_zp                          butter_ap_zp                  ;

p_cheby_poly1                           cheby_poly1                   ;
p_cheby_poly2                           cheby_poly2                   ;
p_cheby1_ap                             cheby1_ap                     ;
p_cheby1_ap_zp                          cheby1_ap_zp                  ;
p_cheby2_ap                             cheby2_ap                     ;
p_cheby2_ap_wp1                         cheby2_ap_wp1                 ;
p_cheby2_ap_zp                          cheby2_ap_zp                  ;
p_cmplx2re                              cmplx2re                      ;
p_concat                                concat                        ;
p_conv                                  conv                          ;
p_conv_cmplx                            conv_cmplx                    ;
p_conv_fft                              conv_fft                      ;
p_conv_fft_cmplx                        conv_fft_cmplx                ;
p_cos_cmplx                             cos_cmplx                     ;

p_decimate                              decimate                      ;
p_decimate_cmplx                        decimate_cmplx                ;
p_dft                                   dft                           ;
p_dft_cmplx                             dft_cmplx                     ;
p_dmod                                  dmod                          ;
p_dspl_info                             dspl_info                     ;

p_ellip_acd                             ellip_acd                     ;
p_ellip_acd_cmplx                       ellip_acd_cmplx               ;
p_ellip_ap                              ellip_ap                      ;
p_ellip_ap_zp                           ellip_ap_zp                   ;
p_ellip_asn                             ellip_asn                     ;
p_ellip_asn_cmplx                       ellip_asn_cmplx               ;
p_ellip_cd                              ellip_cd                      ;
p_ellip_cd_cmplx                        ellip_cd_cmplx                ;
p_ellip_landen                          ellip_landen                  ;
p_ellip_modulareq                       ellip_modulareq               ;
p_ellip_rat                             ellip_rat                     ;
p_ellip_sn                              ellip_sn                      ;
p_ellip_sn_cmplx                        ellip_sn_cmplx                ;

p_farrow_lagrange                       farrow_lagrange               ;
p_farrow_spline                         farrow_spline                 ;
p_fft                                   fft                           ;
p_fft_abs                               fft_abs                       ;
p_fft_abs_cmplx                         fft_abs_cmplx                 ;
p_fft_cmplx                             fft_cmplx                     ;
p_fft_create                            fft_create                    ;
p_fft_free                              fft_free                      ;
p_fft_mag                               fft_mag                       ;
p_fft_mag_cmplx                         fft_mag_cmplx                 ;
p_fft_shift                             fft_shift                     ;
p_fft_shift_cmplx                       fft_shift_cmplx               ;
p_filter_freq_resp                      filter_freq_resp              ;
p_filter_iir                            filter_iir                    ;
p_filter_ws1                            filter_ws1                    ;
p_filter_zp2ab                          filter_zp2ab                  ;
p_find_max_abs                          find_max_abs                  ;
p_find_nearest                          find_nearest                  ;
p_fir_linphase                          fir_linphase                  ;
p_flipip                                flipip                        ;
p_flipip_cmplx                          flipip_cmplx                  ;
p_fourier_integral_cmplx                fourier_integral_cmplx        ;
p_fourier_series_dec                    fourier_series_dec            ;
p_fourier_series_dec_cmplx              fourier_series_dec_cmplx      ;
p_fourier_series_rec                    fourier_series_rec            ;
p_freqs                                 freqs                         ;
p_freqs_cmplx                           freqs_cmplx                   ;
p_freqs2time                            freqs2time                    ;
p_freqz                                 freqz                         ;

p_gnuplot_close                         gnuplot_close                 ;
p_gnuplot_cmd                           gnuplot_cmd                   ;
p_gnuplot_create                        gnuplot_create                ;
p_gnuplot_open                          gnuplot_open                  ;
p_goertzel                              goertzel                      ;
p_goertzel_cmplx                        goertzel_cmplx                ;
p_group_delay                           group_delay                   ;

p_histogram                             histogram                     ;
p_histogram_norm                        histogram_norm                ;

p_idft_cmplx                            idft_cmplx                    ;
p_ifft_cmplx                            ifft_cmplx                    ;
p_iir                                   iir                           ;

p_linspace                              linspace                      ;
p_log_cmplx                             log_cmplx                     ;
p_logspace                              logspace                      ;
p_low2bp                                low2bp                        ;
p_low2bs                                low2bs                        ;
p_low2high                              low2high                      ;
p_low2low                               low2low                       ;

p_matrix_eig_cmplx                      matrix_eig_cmplx              ;
p_matrix_eye                            matrix_eye                    ;
p_matrix_eye_cmplx                      matrix_eye_cmplx              ;
p_matrix_mul                            matrix_mul                    ;
p_matrix_pinv                           matrix_pinv                   ;
p_matrix_print                          matrix_print                  ;
p_matrix_print_cmplx                    matrix_print_cmplx            ;
p_matrix_svd                            matrix_svd                    ;
p_matrix_transpose                      matrix_transpose              ;
p_matrix_transpose_cmplx                matrix_transpose_cmplx        ;
p_matrix_transpose_hermite              matrix_transpose_hermite      ;
p_mean                                  mean                          ;
p_mean_cmplx                            mean_cmplx                    ;
p_minmax                                minmax                        ;

p_ones                                  ones                          ;

p_phase_delay                           phase_delay                   ;
p_poly_z2a_cmplx                        poly_z2a_cmplx                ;
p_polyroots                             polyroots                     ;
p_polyval                               polyval                       ;
p_polyval_cmplx                         polyval_cmplx                 ;
p_psd_bartlett                          psd_bartlett                  ;
p_psd_bartlett_cmplx                    psd_bartlett_cmplx            ;
p_psd_periodogram                       psd_periodogram               ;
p_psd_periodogram_cmplx                 psd_periodogram_cmplx         ;
p_psd_welch                             psd_welch                     ;
p_psd_welch_cmplx                       psd_welch_cmplx               ;

p_randb                                 randb                         ;
p_randb2                                randb2                        ;
p_randi                                 randi                         ;
p_randn                                 randn                         ;
p_randn_cmplx                           randn_cmplx                   ;
p_random_init                           random_init                   ;
p_randu                                 randu                         ;
p_ratcompos                             ratcompos                     ;
p_re2cmplx                              re2cmplx                      ;
p_readbin                               readbin                       ;

p_signal_pimp                           signal_pimp                   ;
p_signal_saw                            signal_saw                    ;
p_sin_cmplx                             sin_cmplx                     ;
p_sinc                                  sinc                          ;
p_sine_int                              sine_int                      ;
p_sqrt_cmplx                            sqrt_cmplx                    ;
p_stat_std                              stat_std                      ;
p_stat_std_cmplx                        stat_std_cmplx                ;
p_sum                                   sum                           ;
p_sum_sqr                               sum_sqr                       ;

p_trapint                               trapint                       ;
p_trapint_cmplx                         trapint_cmplx                 ;

p_unwrap                                unwrap                        ;

p_vector_dot                            vector_dot                    ;
p_verif                                 verif                         ;
p_verif_data_gen                        verif_data_gen                ;
p_verif_cmplx                           verif_cmplx                   ;
p_verif_str                             verif_str                     ;
p_verif_str_cmplx                       verif_str_cmplx               ;

p_window                                window                        ;
p_writebin                              writebin                      ;
p_writetxt                              writetxt                      ;
p_writetxt_3d                           writetxt_3d                   ;
p_writetxt_3dline                       writetxt_3dline               ;
p_writetxt_cmplx                        writetxt_cmplx                ;
p_writetxt_cmplx_im                     writetxt_cmplx_im             ;
p_writetxt_cmplx_re                     writetxt_cmplx_re             ;
p_writetxt_int                          writetxt_int                  ;

p_xcorr                                 xcorr                         ;
p_xcorr_cmplx                           xcorr_cmplx                   ;


#ifdef WIN_OS
#define LOAD_FUNC(fn) \
        fname = #fn;\
        fn = (p_##fn)GetProcAddress(handle, fname);\
        if(! fn) goto exit_label;
#endif



#ifdef LINUX_OS
#define LOAD_FUNC(fn) \
        fname = #fn;\
        fn = (p_##fn)dlsym(handle, fname);\
        if ((error = dlerror()) != NULL) goto exit_label
#endif


void* dspl_load()
{
    char* fname;
    #ifdef WIN_OS
        HINSTANCE handle;
        handle = LoadLibrary(TEXT("libdspl.dll"));
        if (!handle)
        {
            printf("libdspl.dll loading ERROR!\n");
            return NULL;
        }
    #endif /* WIN_OS */
    
    
    #ifdef LINUX_OS
        char* error;
        void *handle;
        /* open the *.so */
        handle = dlopen ("./libdspl.so", RTLD_LAZY);
        if (!handle)
        {
            printf("libdspl.so loading ERROR!\n");
            return NULL;
        }
    #endif  /* LINUX_OS */
    
    LOAD_FUNC(acos_cmplx);
    LOAD_FUNC(addlog);
    LOAD_FUNC(array_scale_lin);
    LOAD_FUNC(asin_cmplx);
    
    LOAD_FUNC(bessel_i0);
    LOAD_FUNC(bilinear);
    LOAD_FUNC(butter_ap);
    LOAD_FUNC(butter_ap_zp);
    
    LOAD_FUNC(cheby_poly1);
    LOAD_FUNC(cheby_poly2);
    LOAD_FUNC(cheby1_ap);
    LOAD_FUNC(cheby1_ap_zp);
    LOAD_FUNC(cheby2_ap);
    LOAD_FUNC(cheby2_ap_wp1);
    LOAD_FUNC(cheby2_ap_zp);
    LOAD_FUNC(cmplx2re);
    LOAD_FUNC(concat);
    LOAD_FUNC(conv);
    LOAD_FUNC(conv_cmplx);
    LOAD_FUNC(conv_fft);
    LOAD_FUNC(conv_fft_cmplx);
    LOAD_FUNC(cos_cmplx);
    
    LOAD_FUNC(decimate);
    LOAD_FUNC(decimate_cmplx);
    LOAD_FUNC(dft);
    LOAD_FUNC(dft_cmplx);
    LOAD_FUNC(dmod);
    LOAD_FUNC(dspl_info);
    
    LOAD_FUNC(ellip_acd);
    LOAD_FUNC(ellip_acd_cmplx);
    LOAD_FUNC(ellip_ap);
    LOAD_FUNC(ellip_ap_zp);
    LOAD_FUNC(ellip_asn);
    LOAD_FUNC(ellip_asn_cmplx);
    LOAD_FUNC(ellip_cd);
    LOAD_FUNC(ellip_cd_cmplx);
    LOAD_FUNC(ellip_landen);
    LOAD_FUNC(ellip_modulareq);
    LOAD_FUNC(ellip_rat);
    LOAD_FUNC(ellip_sn);
    LOAD_FUNC(ellip_sn_cmplx);
    
    LOAD_FUNC(farrow_lagrange);
    LOAD_FUNC(farrow_spline);
    LOAD_FUNC(fft);
    LOAD_FUNC(fft_cmplx);
    LOAD_FUNC(fft_create);
    LOAD_FUNC(fft_free);
    LOAD_FUNC(fft_mag);
    LOAD_FUNC(fft_mag_cmplx);
    LOAD_FUNC(fft_shift);
    LOAD_FUNC(fft_shift_cmplx);
    LOAD_FUNC(filter_freq_resp);
    LOAD_FUNC(filter_iir);
    LOAD_FUNC(filter_ws1);
    LOAD_FUNC(filter_zp2ab);
    LOAD_FUNC(find_max_abs);
    LOAD_FUNC(find_nearest);
    LOAD_FUNC(fir_linphase);
    LOAD_FUNC(flipip);
    LOAD_FUNC(flipip_cmplx);
    LOAD_FUNC(fourier_integral_cmplx);
    LOAD_FUNC(fourier_series_dec);
    LOAD_FUNC(fourier_series_dec_cmplx);
    LOAD_FUNC(fourier_series_rec);
    LOAD_FUNC(freqz);
    LOAD_FUNC(freqs);
    LOAD_FUNC(freqs_cmplx);
    LOAD_FUNC(freqs2time);
    
    LOAD_FUNC(gnuplot_close);   
    LOAD_FUNC(gnuplot_cmd);     
    LOAD_FUNC(gnuplot_create);
    LOAD_FUNC(gnuplot_open);
    LOAD_FUNC(goertzel);
    LOAD_FUNC(goertzel_cmplx);
    LOAD_FUNC(group_delay);
    
    LOAD_FUNC(histogram);
    LOAD_FUNC(histogram_norm);
    
    LOAD_FUNC(idft_cmplx);
    LOAD_FUNC(ifft_cmplx);
    LOAD_FUNC(iir);
    
    LOAD_FUNC(linspace);
    LOAD_FUNC(log_cmplx);
    LOAD_FUNC(logspace);
    LOAD_FUNC(low2bp);
    LOAD_FUNC(low2bs);
    LOAD_FUNC(low2high);
    LOAD_FUNC(low2low);
    
    LOAD_FUNC(matrix_eig_cmplx);
    LOAD_FUNC(matrix_eye);
    LOAD_FUNC(matrix_eye_cmplx);
    LOAD_FUNC(matrix_mul);
    LOAD_FUNC(matrix_pinv);
    LOAD_FUNC(matrix_print);
    LOAD_FUNC(matrix_print_cmplx);
    LOAD_FUNC(matrix_svd);
    LOAD_FUNC(matrix_transpose);
    LOAD_FUNC(matrix_transpose_cmplx);
    LOAD_FUNC(matrix_transpose_hermite);
    LOAD_FUNC(mean);
    LOAD_FUNC(mean_cmplx);
    LOAD_FUNC(minmax);
    
    LOAD_FUNC(ones);
    
    LOAD_FUNC(phase_delay);
    LOAD_FUNC(poly_z2a_cmplx);
    LOAD_FUNC(polyroots);
    LOAD_FUNC(polyval);
    LOAD_FUNC(polyval_cmplx);
    LOAD_FUNC(psd_bartlett);
    LOAD_FUNC(psd_bartlett_cmplx);
    LOAD_FUNC(psd_periodogram);
    LOAD_FUNC(psd_periodogram_cmplx);
    LOAD_FUNC(psd_welch);
    LOAD_FUNC(psd_welch_cmplx);
    
    LOAD_FUNC(randi);
    LOAD_FUNC(randb);
    LOAD_FUNC(randb2);
    LOAD_FUNC(randn);
    LOAD_FUNC(randn_cmplx);
    LOAD_FUNC(random_init);
    LOAD_FUNC(randu);
    LOAD_FUNC(ratcompos);
    LOAD_FUNC(re2cmplx);
    LOAD_FUNC(readbin);
    
    LOAD_FUNC(signal_pimp);
    LOAD_FUNC(signal_saw);
    LOAD_FUNC(sin_cmplx);
    LOAD_FUNC(sinc);
    LOAD_FUNC(sine_int);
    LOAD_FUNC(sqrt_cmplx);
    LOAD_FUNC(stat_std);
    LOAD_FUNC(stat_std_cmplx);
    LOAD_FUNC(sum);
    LOAD_FUNC(sum_sqr);
    
    LOAD_FUNC(trapint);
    LOAD_FUNC(trapint_cmplx);
    
    LOAD_FUNC(unwrap);
    
    LOAD_FUNC(vector_dot);
    LOAD_FUNC(verif);
    LOAD_FUNC(verif_data_gen);
    LOAD_FUNC(verif_cmplx);
    LOAD_FUNC(verif_str);
    LOAD_FUNC(verif_str_cmplx);
    
    LOAD_FUNC(window);
    LOAD_FUNC(writebin);
    LOAD_FUNC(writetxt);
    LOAD_FUNC(writetxt_3d);
    LOAD_FUNC(writetxt_3dline);
    LOAD_FUNC(writetxt_cmplx);
    LOAD_FUNC(writetxt_cmplx_im);
    LOAD_FUNC(writetxt_cmplx_re);
    LOAD_FUNC(writetxt_int);
    
    LOAD_FUNC(xcorr);
    LOAD_FUNC(xcorr_cmplx);

    
    #ifdef WIN_OS
    return (void*)handle;
    exit_label:
        printf("function %s loading ERROR!\n", fname);
        if(handle)
            FreeLibrary(handle);
        return NULL;
    #endif /* WIN_OS */
    
    
    #ifdef LINUX_OS
      return handle;
    exit_label:
        printf("function %s loading ERROR!\n", fname);
        if(handle)
            dlclose(handle);
        return NULL;
    #endif /* LINUX_OS */
}



void dspl_free(void* handle)
{
    #ifdef WIN_OS
        FreeLibrary((HINSTANCE)handle);
    #endif /* WIN_OS */
    
    #ifdef LINUX_OS
        dlclose(handle);
    #endif /* LINUX_OS */
}

#endif /* BUILD_LIB */
