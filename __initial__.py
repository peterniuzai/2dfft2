import numpy as np
import matplotlib.pyplot as plt
import sys,os,time
import readGBT
import pyfits
import h5py
import mpi4py.MPI as MPI
import scipy.signal as signal
from scipy.interpolate import griddata
from readfile import read_data
from dir_create import dir_create
from calibrated import calibration
from rebin import rebin, rebin_inter
from FFT import FFT
from polar_transform_inter import polar_coordinates_convert
from Signal_finding import Signal_finding
from DM_calculate import DM_calculate
from plot_all import plot

import warnings
warnings.filterwarnings("ignore")
Functions to process this method

def distribute_data(d_sets,t_sets,freq,p_n,t_len,datatype):

#Then distribute them to multiprocess, 
#Finally,return the data sets 
     if p_n < 1 :
           print 'the data length is shortter than ',t_len*2,'!'   
           freq = comm.bcast(freq,root=0)
	   d_sets_loc = comm.bcast( d_sets,root =0)
           t_sets_loc = comm.bcast( t_sets,root =0)
     else:  
           freq = comm.bcast(freq,root=0)
           data = np.empty((4096,t_len*p_n),dtype = datatype)
           time = np.empty(t_len*p_n)
           comm.barrier()
	   comm.Scatter(d_sets ,data, root = 0)
#           data = comm.scatter(d_sets , root = 0)
#    data = comm.scatter(d_sets , root = 0)
#           comm.barrier()
#           time = comm.scatter(t_sets , root = 0)
	   time =  comm.scatter(t_sets , root = 0)
           comm.barrier()
           d_sets_loc = []
           t_sets_loc = []
           for i in range(p_n):
                       d_sets_loc.append( data[:, i* t_len : (i+1)* t_len] )
                       t_sets_loc.append( time[   i* t_len : (i+1)* t_len] )


     return d_sets_loc,t_sets_loc,freq


if __name__ == '__main__':

##############################################################################################
#Arguments 
##############################################################################################


     # instance for invoking MPI relatedfunctions
     comm = MPI.COMM_WORLD
     # the node rank in the whole community
     comm_rank = comm.Get_rank()
     # the size of the whole community, i.e.,the total number of working nodes in the MPI cluster
     comm_size = comm.Get_size()

##############################################################################################
# B1929+10_filtered.npy B2319+60_filtered.npy B0329+54_filtered.npy 
# filtered.npy  calibrated.npy  filtered_short.npy high_passed_short.npy  raw.npy
     f_name    = 'B0329+54-25db.fits'
#     f_name    = 'B0329+54-10db.fits'
#     f_name    = 'B1859+03-00db.fits'
#     f_name    = 'B1929+10-00db.fits'
#     f_name  = 'B1859+03-10db.fits'
#     f_name    = 'wiggle.npy'
#     f_name    = str(cut_id) + '.npy'
#     f_name    = str(3) + '.npy'
#     f_name    = '3_1.npy'
#     f_name    = 'B0329+54_filtered.npy'
#     f_name    = 'filtered_short.npy'
#     f_name    = 'filtered.npy'
#     f_name    =  'wiggle_J2139+00h5'
#     f_dir     =  '/home/nch/pulsar_data/wiggle_cut/'
     f_dir     = '/home/nch/Data/pulsar_data/'
#     f_dir	= '/home/nch/Data/burst_data/'
#     f_dir     = '/home/nch/FFT_search/data/wigglez_found/data_pick/'
####################################################################################################

     t_len     = 2048			#time length for each smallest unit to process.
     rad_grid  = 1024			#radius grid size for interpolate in polar-coordin transform.
     ang_grid  = 860			#angle grid size for interpolate in polar-coordin transform.
     ang_min   = -85			#range of angle in polar transform :minum value.
     ang_max   = 0			#range of angle in polar transform :max value.
     msk_cycle = 5			#the number of channels to be zeros in 2D-FFT(Noise remove).
     pixel     = 2			#the number of pixel to sum in 2ndFFT3D SNR compute.

 
##############################################################################################
#start to process
##############################################################################################
     time_1    = time.time()

     if comm_rank == 0: 
		print 'Begin to load data from ' + f_name + '....:',comm_rank
		d_sets , t_sets, freq, p_n ,datatype,plot_file = read_data(f_dir, f_name , t_len, comm_size)
     else:
	  	datatype = p_n = freq = d_sets = t_sets = plot_file= None

#####################################################################################################
#plot arguments
##################     
     plot_file   = comm.bcast( plot_file,root=0)
     p_n         = comm.bcast( p_n,root=0)
     datatype    = comm.bcast( datatype,root=0)
     comm.barrier()
     plot_dir  = '../graph/' + plot_file + '/'
     plot_proc = 'raw,rebin,1stFFT,2ndFFT_2D,2ndFFT_3D,polar_sets_3D,polar_sets_2D'
 #   plot_proc = 'raw'
     if comm_rank == 0:
                dir_create(plot_dir)
                print 'directory build complete!'
##############################################################################################
     if comm_rank == 0: print 'load data over ,begin to distribute them to multiprocess!'               
     d_sets, t_sets, freq = distribute_data(d_sets, t_sets,freq ,p_n, t_len, datatype)
     if comm_rank == 0:    print 'Data load and scatter to multiprocess over. rank:',comm_rank

#     if comm_rank == 0:    print 'Begin to calibrate... rank:',comm_rank
#     chan_equaliz = np.load('/home/nch/FFT_search/src/chan_equaliz.npy')
#     d_sets = calibration(chan_equaliz,d_sets)

     if comm_rank == 0:    print 'Begin to rebin... rank:',comm_rank
     re_sets, f_axis ,nbin =  rebin(d_sets,t_sets, freq)

     if comm_rank == 0:    print 'Rebin over. rank:',comm_rank

     if comm_rank == 0:    print 'Begin to do 1st 2-D FFT on rebin data... rank :',comm_rank
     FFT1st_sets = FFT(re_sets, 2,msk_cycle,comm_rank)

     if comm_rank == 0:    print '1st FFT over. rank :',comm_rank

     if comm_rank == 0:    print 'Begin to transform rectangular coordinates into polar coordinates.....  rank:',comm_rank
     polar_sets  = polar_coordinates_convert( FFT1st_sets,rad_grid ,ang_grid, ang_min ,ang_max )

     if comm_rank == 0:    print 'Polar transform over. rank:',comm_rank

     if comm_rank == 0:    print 'Begin to do the 2nd 1-D FFT along radius direction.... rank:',comm_rank
     FFT2nd_sets = FFT( polar_sets, 1 )

     if comm_rank == 0:    print '2nd FFT over. rank:',comm_rank

     if comm_rank == 0:    print 'Begin to locate the signal and calculate SNR... rank:',comm_rank
     SNR , location = Signal_finding(FFT2nd_sets,  ang_min , ang_max, pixel)

     if comm_rank == 0:    print 'Calculation over. rank:',comm_rank
     time_2  =  time.time()
     consume =  time_2 - time_1

     if comm_rank == 0:    print str(comm_rank)+' process matrix size: 4096 * '+str(len(SNR)*2048) + '\ntime is:', consume ,'seconds,  ','equal',consume/60.,'minutes.'

     if comm_rank == 0:    print 'start to plot... rank:',comm_rank
     plot(comm_rank,t_sets,d_sets,re_sets,polar_sets,FFT1st_sets,FFT2nd_sets,plot_proc,freq,f_axis,2,rad_grid,ang_grid,plot_dir,pixel,ang_min,ang_max)

###############################################################################################
#gather the results from all processes
##############################################################################################
     comm.barrier()
     combine_SNR      = comm.gather(SNR,root=0)
     combine_SNR      = np.array(combine_SNR).reshape(-1)
     combine_location = comm.gather(location,root=0)
     combine_location = np.array(combine_location).reshape(-1)

     if comm_rank == 0:   
                print 'multiprocess plot over....'
                N_cut = len(combine_SNR)
                plt.ylabel('SNR')
                plt.xlabel(str(N_cut)+'files')
                plt.title('SNR of each '+str(N_cut)+' files')
                plt.plot(np.arange(N_cut),combine_SNR,label='SNR of 2nd FFT in 2-D map',color='blue')
                plt.plot(np.arange(N_cut),combine_SNR,'ro')
#               plt.legend(loc='upper left')
                plt.grid()
                p_dir = plot_dir + 'SNR/'
                plt.savefig(plot_dir+'SNR_'+ str(rad_grid)+'_'+str(ang_grid))
#                plt.show()
                plt.close()

                N_cut = len(combine_location)
                plt.ylabel('Location(degree)')
                plt.xlabel(str(N_cut)+'files')
                plt.title('Location of each '+str(N_cut)+' files')
                plt.plot(np.arange(N_cut),combine_location,label='Location of 2nd FFT in 2-D map',color='blue')
                plt.plot(np.arange(N_cut),combine_location,'ro')
#                plt.legend(loc='upper left')
                plt.grid()
                p_dir = plot_dir + 'Location/'
                plt.savefig(plot_dir + 'Loc_'+ str(rad_grid)+'_'+str(ang_grid))
#                plt.show()
                plt.close()

       	        print 'plot over.... rank:',comm_rank

