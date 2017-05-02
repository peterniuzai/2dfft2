import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal
def plot(comm_rank,t_sets , d_sets , re_sets , polar_sets , FFT1st_sets , FFT2nd_sets,\
         process , freq , f_axis ,  Dim , r_g , a_g ,dir , pixel = 5, ang_min = -80 , ang_max = 0):



       N_cut = len(d_sets)# Nomber of 2048*2048 tile in each process.

       if   process  == '':
            if comm_rank == 0:    print 'Please claim which process to plot!'
       if 'raw' in process:
            for i in range(N_cut):
                t_x   = t_sets[i]
                data  = np.ma.masked_invalid(np.abs(d_sets[i]))
                mean  = np.mean(data)
                sigma = np.var(data)
                max   = mean + 2 * sigma
                min   = mean - 2 * sigma
                seq   = comm_rank * N_cut + i
                plt.pcolormesh(t_x,freq,data)#,vmax = max,vmin = min)
                plt.title('raw_data')
                plt.xlabel('time(s)')
                plt.ylabel('frequency(Mhz)')
                plt.xlim(t_x.min(),t_x.max())
                plt.ylim(freq.min(),freq.max())
                plt.colorbar()
                p_dir = dir + 'raw/'
                plt.savefig(p_dir + 'raw_'+str(seq))
                plt.close()
       if comm_rank == 0: print 'raw data plot is over...'

       if 'rebin' in process:
            for i in range(N_cut):
                t_x   = t_sets[i]
                data  = np.abs(re_sets[i])
                mean  = np.mean(data)
                sigma = np.var(data)
                max   = mean + 2 * sigma
                min   = mean - 2 * sigma
                seq   = comm_rank * N_cut + i
                plt.pcolormesh(t_x,f_axis,data)#,vmax = max,vmin = min)
                plt.title('data after rebin')
                plt.xlabel('time(s)')
                plt.ylabel('frequency(Mhz)')
                plt.xlim(t_x.min(),t_x.max())
                plt.ylim(f_axis.min(),f_axis.max())
                plt.colorbar()
                p_dir = dir + 'rebin/'
                plt.savefig(p_dir + 'rebin_' + str(seq))
                plt.close()
       if comm_rank==0: print 'rebin plot over...'

       if '1stFFT' in process:
            for i in range(N_cut):
                data  = np.abs(FFT1st_sets[i])
                mean  = np.mean(data)
                sigma = np.var(data)
                max   = mean + 2 * sigma
                min   = mean - 2 * sigma
                seq   = comm_rank * N_cut + i
                plt.pcolormesh(data)#,vmax = max,vmin = min)
                plt.title('1st FFT')
                plt.xlabel('Time axis after 1st FFT')
                plt.ylabel('Frequency after 2nd FFT')
                plt.colorbar()
                index  = np.where(data == np.max(data))
                cord =(index[1][0],index[0][0])
#                if comm_rank == 0:    print index
#                np.savetxt('/home/nch/cord.txt',cord)
                plt.annotate('max:'+str(index[1][0]), xy = cord, xytext = cord)
                              #arrowprops = dict(facecolor = 'red', shrink = 0.01))
                p_dir = dir + '1stFFT/'
                plt.savefig(p_dir + '1stFFT_' + str(seq))
                plt.close()
       if comm_rank==0: print '1stFFT plot over....'

       if 'polar_sets_3D' in process:
            for i in range(N_cut):
                data  = np.abs(polar_sets[i])
                mean  = np.mean(data)
                sigma = np.var(data)
                max   = mean + 2 * sigma
                min   = mean - 2 * sigma
                seq   = comm_rank * N_cut + i
                SNR   = (data.max()-data.mean())/data.std()
                x_axis  = np.linspace(ang_min,ang_max,data.shape[1])
                y_axis  = np.arange(data.shape[0])
                plt.pcolormesh(x_axis,y_axis,data)#,vmax = max,vmin = min)
                plt.title('radius - angle(grid size:'+ str(r_g)+'*'+str(a_g)+')')
                plt.xlabel('Angle(in degree)')
                plt.ylabel('Radius')
                plt.figtext(0.08,0.98,'SNR:'+str(SNR))
                plt.xlim(x_axis.min(),x_axis.max())
                plt.ylim(y_axis.min(),y_axis.max())
                plt.colorbar()
                p_dir = dir + 'polar_sets_3D/'
                plt.savefig(p_dir + 'polar_3D_' + str(seq) )
                plt.close() 
       if comm_rank==0: print 'polar_3D plot over...'

       if 'polar_sets_2D' in process:
            for i in range(N_cut):
                data  = np.abs(polar_sets[i])
                data  = data.sum(axis = 0)
                seq   = comm_rank * N_cut + i
#                if seq == 7:
#                    np.save('/home/nch/temraw.npy',data)
                x_axis  = np.linspace(ang_min,ang_max,data.shape[0])
                #Filter the profile of the FRB signal
                prof_data = signal.medfilt(data,199)
                data      = data - prof_data
                SNR       = (data.max()-data.mean())/data.std()
                dmax =  np.argmax(data)
                cord =  (x_axis[dmax] ,data[dmax])
                plt.plot([cord[0]],[data[dmax]],'ro')
                plt.figtext(0.08,0.98,'SNR:'+str(SNR))
                plt.title('polar Sum along radius axis(grid size:'+ str(r_g)+'*'+str(a_g)+')')
                plt.xlabel('Angle(in degree)')
                plt.ylabel('Intensity')


                plt.annotate('angle:'+str(cord[0])+'deg', xy = cord, xytext = cord, \
                              arrowprops = dict(facecolor = 'black', shrink = 0.1))

                plt.grid()
                plt.xlim(x_axis.min(),x_axis.max())
                plt.ylim(data.min()-10,data.max()+10)
                plt.plot(x_axis,data)
                p_dir = dir + 'polar_sets_2D/'
                plt.savefig(p_dir + 'polar_2D_' + str(seq))
                plt.close()
       if comm_rank==0: print 'polar_2D plot over...'

       if '2ndFFT_3D' in process:
            for i in range(N_cut):
                data  = np.abs(FFT2nd_sets[i])
                mean  = np.mean(data)
                sigma = np.var(data)
                max   = mean + 2 * sigma
                min   = mean - 2 * sigma
                seq   = comm_rank * N_cut + i
                lo    = np.where(data == np.max(data))
                d_max = 0
                for i in np.arange(-pixel,pixel):
                    for j in np.arange(-pixel,pixel):
                        d_max += data[lo[0][0]+i,lo[1][0]+j]
                SNR   = (d_max - data.mean())/data.std()
                for ii in range(3):
                       ind  = np.where(data == data.max())
                       y_max  = ind[0][0]
                       if y_max  ==  data.shape[0]/2 or y_max == data.shape[0]/2 + 1:
                          data[ind] = 0
                ind  = np.where(data == data.max())
                y_ax  = (ind[0][0]-data.shape[0]/2)
                x_axis  = np.linspace(ang_min,ang_max,data.shape[1])
                deg   = x_axis[(ind[1][0])]

                y_axis  = np.arange(-data.shape[0]/2,data.shape[0]/2)
                plt.pcolormesh(x_axis,y_axis,data)#,vmax = max,vmin = min)
                plt.title('2nd FFT along radius axis(grid size:'+ str(r_g)+'*'+str(a_g)+')')
                plt.xlabel('Angle(in degree)')
                plt.ylabel('Radius after FFT')
                plt.xlim(x_axis.min(),x_axis.max())
                plt.ylim(y_axis.min(),y_axis.max())
                plt.figtext(0.08,0.98,'SNR:'+ str(int(SNR)) + ' Location:' + str(deg) + 'y_axis'+str(y_ax))
                plt.colorbar()
                p_dir = dir + '2ndFFT_3D/'
                plt.savefig(p_dir +'2ndFFT_3D_' + str(seq) )
                plt.close()
       if comm_rank==0: print '2nd FFT 3D plot over...'

       if '2ndFFT_2D' in process:
            for i in range(N_cut):
                data  = np.abs(FFT2nd_sets[i])
                data  = data.sum(axis = 0)
                seq   = comm_rank * N_cut + i
                x_axis    = np.linspace(ang_min,ang_max,data.shape[0])
                #Filter the profile of the FRB signal
                prof_data = signal.medfilt(data,199)
                data      = data - prof_data
                SNR       = (data.max()-data.mean())/data.std()
                dmax =  np.argmax(data)
                cord =  (x_axis[dmax] ,data[dmax])
                plt.plot([cord[0]],[data[dmax]],'ro')
                plt.title('Sum along radius axis(grid size:'+ str(r_g)+'*'+str(a_g)+')')
                plt.xlabel('Angle(in degree)')
                plt.ylabel('Intensity')

                plt.annotate('angle:'+str(cord[0])+'deg', xy = cord, xytext = cord, \
                              arrowprops = dict(facecolor = 'black', shrink = 0.1))

                plt.grid()
                plt.xlim(x_axis.min(),x_axis.max())
                plt.ylim(data.min()-10,data.max()+10)
                plt.figtext(0.08,0.98,'SNR:'+str(SNR))
                plt.plot(x_axis,data)
                p_dir = dir + '2ndFFT_2D/'
                plt.savefig(p_dir + '2ndFFT_2D_' + str(seq))
                plt.close()
       if comm_rank==0: print '2ndFFT 2D plot over...'

if __name__ == '__main__':
                exit()          



