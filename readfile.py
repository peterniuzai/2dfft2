import numpy as np
import readGBT
import pyfits
import h5py
def read_data(f_dir, f_name, t_len = 2048,comm_size=10):

#load data from GBT, chunked them to data sets with shape :t_length * t_length (Default :2048*2048)
#Then distribute them to multiprocess, 
#Finally,return the data sets 
         if '.fits' in f_name :
         		hdulist = pyfits.open(f_dir + f_name)
		        data ,tx, ra, dec, az, el, freq = readGBT.read_fits(hdulist)
                        if data.shape[1] == 4:
		                 data   = data[:,0,:]
		        else:
                		 data   = data
			plot_file = f_name[:-5]

                     
         if 'h5'    in f_name : 
           		file   = h5py.File(f_dir + f_name ,'r')
		        data   = file['DATA'][:]
                        out = np.empty((4096, 4, 1003 * 2048), dtype=data.dtype)
		        for ii in range(1003):
			        for jj in range(4):
			                 out[:,jj,2048 * ii:2048 * (ii + 1)] = data[ii,:,jj,:].T
			data = out
                        if data.shape[1] == 4:
		                 data   = data[:,0,:]
		        else:
                		 data   = data

			tx     = np.arange(data.shape[1])
			plot_file = f_name[:-2]
	
         if '.npy'  in f_name:
			data   = np.load(f_dir + f_name)
                        if data.shape[1] == 4:
		                 data   = data[:,0,:]
				 print 'data shape is [:,4,:]'
		        else:
                		 data   = data
				 print 'data shape is [:,0,:]'
    		        tx     = np.arange(data.shape[1])    
			plot_file = f_name[:-4]

         freq   = np.load(f_dir + 'freq.npy')

         num    = data.shape[1]/t_len
         p_n    = num/comm_size
         d_sets = []
         t_sets = []
	 datatype =data.dtype
	 data = np.nan_to_num(data)
         if p_n < 1:
		   print 'short data~'
         	   d_sets = [data]
                   t_sets = [tx]
         else:
	           for i in range(comm_size):
          	        	 d_sets.append( data[:, i* t_len * p_n : (i+1)* t_len * p_n] )
  	                	 t_sets.append(   tx[   i* t_len * p_n : (i+1)* t_len * p_n] )
         d_sets = np.array(d_sets)
         t_sets = np.array(t_sets)

         return d_sets, t_sets, freq, p_n, datatype, plot_file

if   __name__ == '__main__':
     f_dir  = '/home/nch/pulsar_data/'
     f_name = 'B1929+10-00db.fits'
     d ,t ,freq ,p_n = read_data(f_dir ,f_name )
     print 'load over!'
     exit()
