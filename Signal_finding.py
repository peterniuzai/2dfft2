import numpy as np
import matplotlib.pyplot as plt
import sys,os,time
import scipy.signal as signal
def Signal_finding(FFT2_sets, ang_min = -90,ang_max = 0, pixel=2):
    SNR        = []
    location   = []
    for i in range(len(FFT2_sets)):
         data  = FFT2_sets[i]
         deg   = np.linspace(ang_min,ang_max,data.shape[1])
         data  = abs(data)
         lo    = np.where(data == np.max(data))
         d_max = data.max()
#         for i in np.arange(-pixel,pixel):
#             for j in np.arange(-pixel,pixel):
#                d_max += data[lo[0][0]+i,lo[1][0]+j]
         snr   = (d_max - data.mean())/data.std()

         indx  = lo[1][0]
         SNR.append(snr)
         location.append(deg[indx])
    SNR = np.array(SNR)
    location = np.array(location)
    return SNR , location

if __name__ == '__main__':
        ang_min = -90
        ang_max = 0
        pixel = 3
        data  = np.load('/home/ycli/data/burst_data/filtered_short.npy')
        data  = np.load('/home/nch/plot_raw/wigglez_found/data_pick/0_0.npy')
        if   data.shape[1] == 4:
             data = [data[:,0,:]]
        else:
             data = [data]
        SNR , location = Signal_finding(data, pixel)
        print SNR,location
