import numpy as np
import matplotlib.pyplot as plt
import sys,os,time
import scipy.signal as signal
#import readGBT
#import pyfits
#import h5py
#import mpi4py.MPI as MPI

import scipy.signal as signal
from scipy.interpolate import griddata
#from readfile import read_data
#from dir_create import dir_create
#from calibrated import calibration
#from rebin import rebin, rebin_inter
#from FFT import FFT
#from polar_transform_inter import polar_coordinates_convert
#from Signal_finding import Signal_finding
#from DM_calculate import DM_calculate
#from plot_all import plot
#f_n     = 'FRB010125'
f_n      = 'signal_sim_snr_1_c'
#f_n     = 'FRB110220'
#f_n     = 'FRB010621_5'
#f_n     = sys.argv[1]
#f_q     = 'FRB010125'
#f_q     = 'FRB110220'
#f_q     = 'FRB010621'

#f_q     = 'freq'
f_dir   = '../data/'
#data = np.load(f_dir + 'Crab_52277_022_short.npy')
#data    = np.load(f_dir + 'signal_sim_snr1.npy')
#data    = np.load(f_dir + 'FRB110220.npy')
#data    = np.load(f_dir + 'FRB010125.npy')
data    = np.load(f_dir + f_n + '.npy')
#print data.shape,'2222'
#exit()
#data[1000:1200,:]=np.random.random((200,data.shape[1]))
#data[-200:,:]= np.random.random((200,data.shape[1]))
#data[:200,:]= np.random.random((200,data.shape[1]))
#data    = data[:,:256]
#T = 2
#F = 3000/T
#Value = 30
#data[1800:2400,:] = Value*np.random.random((600,256))

#data[1000:1003,170:170+100] = Value*np.random.random((3,100))

#data[100:100+F,17:17+T] = Value*np.random.random((F,T))

#data[3000:3010,200:] = Value*np.random.random((10,100))
#data[1500:1520,80:140] = Value*np.random.random((20,60))
#data[0:20,100:160] = Value*np.random.random((20,60))
#data[1200:1220,200:260] = Value*np.random.random((20,60))
#data[3000:3650,:] = 100000
#data[:700,:] = 100000
#data[:,30:70] = 10000000
#data[:,800:900] = 100000000
#data    = np.random.random(data.shape)*10
pixel = 1
area  = (pixel*2)**2
#print data.shape,data.max()
freq    = np.load(f_dir +'freq_s.npy')
#freq   = np.load(f_dir + 'Crab_52277_022_freq.npy')
#freq    = np.load(f_dir + 'freq.npy')
#freq    = np.load(f_dir + 'FRB110220_freq.npy')
#freq    = np.load(f_dir + 'FRB010125_freq.npy')
#f_n     = 'FRB010125'
fy       = freq**-2
#nbin        = max(data.shape)
nbin        = data.shape[0]
BW	= freq.max() - freq.min()
t_rsl	= 1 #ms
T	= data.shape[1]*t_rsl #ms
C	= 4.15e6 #cm^3 * pc^-1 * MHz**2 * ms
#nbin	= 512
#print data.shape,nbin
#nbin         = 1024
f_axis       = np.linspace(fy.min(),fy.max(),nbin)
b_aray,b_edg = np.histogram(fy,bins=nbin)

data0 = data 
data_0nan = np.nan_to_num(data)
data1 = np.zeros((nbin,data.shape[1]),data.dtype)
for i in np.arange(nbin):
    length  =  b_aray[i]
    if  i == 0:
        index = 0
    elif i == 1:
        index = b_aray[0]
    else:
        index = b_aray[0:i].sum()
    for ii in np.arange(length):
        data1[i,:]+=data_0nan[ii+index,:]
tem = np.ones(data.shape[1])
bin_weight = b_aray[:,None] * tem[None,:]
data1 = data1 / 1.0 / bin_weight
data1 = np.nan_to_num(data1)
data  = data1

#print data.shape
fft2   = np.fft.fft2(data)
fft2   = np.nan_to_num(fft2)
shift  = np.fft.fftshift(fft2)/np.sqrt(fft2.size)
data   = shift[1:shift.shape[0]/2+1,shift.shape[1]/2:]
#data   =  shift[:shift.shape[0]/2,shift.shape[1]/2:]
#print data.shape
#exit()
data[-1:, :] = 0
data[  :,:1] = 0
#data   = data[:,1:]


#print data.shape,'/*/*/*/*/*'
data2  = data

rang   = data.shape
row    = np.arange(-rang[0]+1,1,dtype=np.float32)
rank   = np.arange(rang[1],dtype=np.float32)
print row.max(),row.min(),rank.max(),rank.min()
#exit()
#angle  = np.nan_to_num(np.arctan(row[:,None]/rank[None,:]))#/ np.pi * 180
#angle  = np.tan(angle)
tangent  = np.nan_to_num((row[:,None]/1.0/rank[None,:]))
tangent[:,0] = tangent[:,1]
print tangent.max(),tangent.min(),'*****',tangent.shape
#angle	= np.tan(angle)
#print angle[:,0]
#exit()

radius = np.sqrt(row[:,None]**2 + rank[None,:]**2)

tan    = tangent.reshape(-1)
rad    = radius.reshape(-1)

#ang_max = np.tan(-5/180.*np.pi)
#ang_min = np.tan(-85/180.*np.pi)
ang_max = -2
ang_min = -89
tan_max = np.tan(ang_max/180.*np.pi)
tan_min = np.tan(ang_min/180.*np.pi)
#print tan_max,tan_min 
#exit()
#rad_grid = 512#data.shape[0]*1
#ang_grid = 512# ((ang_max-ang_min)/0.224) #data.shape[1]*1#
#print ang_grid,'$$$$$$$$$$$$$'
points = (rad,tan)
data   = data.reshape(-1)
#grid_r, grid_a = np.mgrid[3:rang[0]:rad_grid*1j, ang_min:ang_max:ang_grid*1j]
#print grid_r.shape, grid_a.shape,'//////////////'
#print grid_a[0,:]
#print grid_a[-1,:]
#exit()
row_rsl	 = (fy.max() - fy.min())**-1
ran_rsl	 = T**-1
DM_rsl	 = 1.205e-7 * t_rsl * freq.max() ** 3 / BW
print DM_rsl

k_rsl	 = DM_rsl * (ran_rsl / row_rsl)* C * nbin
print k_rsl
tan_grid = abs((tan_max-tan_min)/ k_rsl)
rad_grid = data.size / tan_grid
print tan_grid,'tan_grid number',rad_grid
#exit()
grid_r = np.linspace(0,rang[0],int(rad_grid))
#grid_a = np.linspace(ang_min,ang_max,int(ang_grid))
#grid_t = np.tan(grid_a/180.*np.pi)
grid_t	= np.linspace(tan_min,tan_max,int(tan_grid))
#grid_a = np.tan(grid_a/180.*np.pi)
#print grid_r.max(),grid_r.min(),'/*/*/*'
grid_t,grid_r = np.meshgrid(grid_t,grid_r)
#print grid_r[0,:]
#print grid_r[-1,:]
#print grid_r.shape,'88888888'



polar_matrix   = griddata(points,data,(grid_r,grid_t),method='nearest')
polar_matrix   = np.nan_to_num(polar_matrix)
#if rang[0] > rang[1]:
#    short =  rang[1]
#else:
#    short = rang[0]
#polar_matrix[:,-3:] = 0
#polar_matrix[:10,:] = 0
#polar_matrix = polar_matrix[:short,:]
data3 = polar_matrix
##temp:
#np.save('/home/nch/work/ipython_exercise/fft',data2)
#np.save('/home/nch/work/ipython_exercise/ang',ang)
#np.save('/home/nch/work/ipython_exercise/rad',rad)

#print data3.shape,'***data3.shape'
data4 = abs(data3).sum(axis=0)
#data4[-10:] = data4[:10] = data4[10:-5].mean()

#fit
prof_data  = signal.medfilt(data4,111)
data4      = data4 - prof_data


fft1 = np.fft.fft(data3,axis=0)
fft1 = fft1/np.sqrt(fft1.shape[0])
shift  = np.fft.fftshift(fft1,axes=0)
data5 = shift

snr  = (data4.max()-data4.mean())/data4[:50].std()
print 'data4 SRN:',snr,data4.std()
data5=abs(data5)
c=data5
s = ['rx','r+','ro','r>','r<','rv','r^','rs','r.','r*','bx','b+','bo','b>','b<','bv','b^','sb','b.','b*','gx','g+','go','g>','g<','gv','g^','sg','g.','g*']
temx=np.linspace(ang_min,ang_max,c.shape[1])
temy=np.linspace(-c.shape[0]/2,c.shape[0]/2,c.shape[0])
d_maxx=0
for i in range(area*2):
    lo    = np.where(data5 == np.max(data5))
    print 'i:',i,lo,'Value:',data5[lo[0][0],lo[1][0]]
    d_maxx = d_maxx + np.max(data5)
    data5[lo[0][0],lo[1][0]]=0
#    plt.plot(temx[lo[1][0]],temy[lo[0][0]],s[i],label=str(i)+':'+str(data5.max()))

#plt.legend(loc='best')

#plt.show()
#print c.max(),c.mean(),c.std()
#exit()
#data5 = c
#lo    = np.where(data5 == np.max(data5))
#d_max = 0
SNRt = (data5.max()-data5.mean())/data5.std()
#for i in np.arange(-pixel,pixel):
#    for j in np.arange(-pixel,pixel):
#        d_max += data5[lo[0][0]+i,lo[1][0]+j]
SNR  = (d_maxx-data5[:,:50].mean())/((data5[:,:50].std())*area)



plt.figure(figsize=(21,14))
plt.subplot(3,2,1)
plt.title('raw data')
plt.xlabel('time')
plt.ylabel('Frequency')
x=np.arange(data0.shape[1])
plt.pcolormesh(x,freq,abs(data0))
plt.colorbar()

plt.subplot(3,2,2)
plt.title('rebin data')
plt.xlabel('time')
plt.ylabel('wave^2')
plt.pcolormesh(abs(data1))
plt.colorbar()

plt.subplot(3,2,3)
plt.title('2D-FFT ')
plt.xlabel('time axis after FFT')
plt.ylabel('wave^2 axis after FFT')
plt.pcolormesh(abs(data2))
plt.colorbar()

plt.subplot(3,2,4)
plt.title('Polar coordinates transform')
plt.xlabel('angle')
plt.ylabel('radius')
x_axis = np.linspace(tan_min,tan_max,data3.shape[1])
#x_axis = np.linspace(-85,-5,data3.shape[1])
#x_axis = np.tan(x_axis/180*np.pi)
#x_axis = np.tan(x_axis/180.*np.pi)
y_axis = np.arange(data3.shape[0])
plt.pcolormesh(x_axis,y_axis,abs(data3))
plt.colorbar()

plt.subplot(3,2,5)
plt.title('sum along radius at each angle(SNR:'+str(int(snr))+')')
plt.xlabel('angle from '+str(ang_min)+' to '+str(ang_max)+' degree')
plt.ylabel('abs value of sum at each angle')
plt.plot(x_axis,data4,'b',label='1st FFT')
p=abs(data5).sum(axis=0)
plt.plot(x_axis,p,'r',label='2nd FFT')
plt.legend(loc='best')

plt.subplot(3,2,6)
x_axis = np.linspace(ang_min,ang_max,data5.shape[1])
y_axis = np.arange(-data5.shape[0]/2,data5.shape[0]/2)
plt.xlabel('SNR:'+str(SNR))
#p=abs(data5).sum(axis=0)
#plt.plot(x_axis,p)
plt.pcolormesh(x_axis,y_axis,abs(data5))
plt.colorbar()
plt.savefig(f_n)
plt.show()
#print 'F:',F,"and  T:",T, '  Value:',Value
#print 'pixel:',F*T
#print 'snr:',snr,'  SNR:',SNR,' SNRt',SNRt
