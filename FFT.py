import numpy as np


def FFT(d_sets, Dim = 2 , msk_cycle = 5 ,comm_rank=0):
    if Dim == 2:   #1st 2-D FFT
       FFT_sets   = []
       for i in range(len(d_sets)):
           data   = d_sets[i]
           fft2   = np.fft.fft2(data)
           fft2   = np.nan_to_num(fft2)

           shift  = np.fft.fftshift(fft2)/np.sqrt(fft2.size)

#           data   = shift[:shift.shape[0]/2-1,2:]
           data   =  shift[:shift.shape[0]/2,shift.shape[1]/2:]#oringinal
           data[-2:, :] = 0
           data[ :,0:2] = 0

           for i in np.arange(msk_cycle):
               x_sum   =  np.abs(data[-100:,:]).sum(axis = 0)
               y_sum   =  np.abs(data[    :,:]).sum(axis = 1)
               x_max   =  np.argmax(x_sum)
               y_max   =  np.argmax(y_sum)
               data[:,	  x_max] = 0
               data[y_max,  :  ] = 0
#               y_max      =  np.argmax(data_short)/data_short.shape[1]
#               x_max      =  np.argmax(data_short)%data_short.shape[1]
#               if data.shape[0]-20  < y_max:
#                       data[-40:,x_max-1:x_max+2]=0
#               if comm_rank == 0:    print '***y_max',y_max,'*****x_max:' ,x_max,' ** data.shape:',data.shape
#           data = data[512:,:512]
#           if comm_rank == 0:    print data.shape,'**2***\n****'
           FFT_sets.append(data)
    elif Dim == 1: #2nd_1-D_FFT along radius
       FFT_sets   = []
       for i in range(len(d_sets)):
           data   = d_sets[i]
           fft1   = np.fft.fft(data,axis=0)
           fft1   = fft1/np.sqrt(fft1.shape[0])
           shift  = np.fft.fftshift(fft1,axes=0)
           FFT_sets.append(shift)
    return FFT_sets

if __name__ == '__main__':
      
        msk_cycle = 2
        data  = np.load('/home/ycli/data/burst_data/filtered_short.npy')
#        data  = np.load('/home/nch/plot_raw/wigglez_found/data_pick/0_0.npy')
        if   data.shape[1] == 4:
             data = [data[:,0,:]]
        else:
             data = [data]
        re_sets  = FFT(data, 2 , msk_cycle)
        print re_sets[0].shape 

