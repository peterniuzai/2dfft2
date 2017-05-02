import scipy.signal as signal
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")
def polar_coordinates_convert(FFT_sets,rad_grid = 2400,ang_grid = 600 , ang_min = -90 ,ang_max = 0 ):
    #Transform the rectanular coordinates into polar coordinates

    polar_sets = []
    for i in range(len(FFT_sets)):
          rang  = FFT_sets[i].shape
          row    = np.arange(-rang[0]+1,1,dtype=np.float32)  # move the center to middle of the matrix
          line   = np.arange(rang[1],dtype=np.float32)
          angle  = np.nan_to_num(np.arctan(row[:,None]/line[None,:]) / np.pi * 180) # calculate the polar coordinates  matrix 
          radius = np.sqrt(row[:,None]**2 + line[None,:]**2)
          ang    = angle.reshape(-1)
          rad    = radius.reshape(-1)

          points = (rad,ang)
          data   = FFT_sets[i].reshape(-1)
#          grid_r, grid_a = np.mgrid[0:rad.max():rad_grid*1j, -85:0:ang_grid*1j]
          grid_r, grid_a = np.mgrid[0:rad.max():rad_grid*1j, ang_min:ang_max:ang_grid*1j]
          polar_matrix   = griddata(points,data,(grid_r,grid_a),method='nearest')
          polar_matrix   = np.nan_to_num(polar_matrix)
          if rang[0] > rang[1]:
               short =  rang[1]
          else:
               short = rang[0]
	  polar_matrix = polar_matrix[:short,:]
          polar_matrix[:,-3:] = 0
          polar_matrix[:10,:] = 0
#          polar_matrix = np.clip(abs(polar_matrix),polar_matrix.min(),24)
          polar_sets.append(polar_matrix)
    return  polar_sets


if __name__ == '__main__':

        rad_grid = 2400
        ang_grid = 600
        data  = np.load('/home/ycli/data/burst_data/filtered_short.npy')
        data  = np.load('/home/nch/plot_raw/wigglez_found/data_pick/0_0.npy')
        if   data.shape[1] == 4:
             data = [data[:,0,:]]
        else:
             data = [data]
        polar_sets  = polar_coordinates_convert( data,rad_grid ,ang_grid )
    
        print polar_sets[0].shape 


