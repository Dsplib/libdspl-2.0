# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 14:48:55 2020

@author: bakhu
"""

import numpy as np
import time
import scipy.fftpack as fft

N =  4194304

x = np.random.randn(N) + 1j * np.random.randn(N)


k = 4
s = N
mflops = np.zeros(22)
for q in range(22):
    z = np.copy(x[0:s])
    t = time.time()
    for i in range(k):
        y = fft.fft(z)
  
    dt = (time.time() - t)/float(k)
    mflops[21-q] = 5.0 * s * np.log2(s) / dt / 1E6 

    print("size = %8d   Mflops = %12.3f" % (s, mflops[21-q]))
    k = int(k * 1.7)
    s = int(s / 2)
