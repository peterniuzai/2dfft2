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
#f_n     = 'FRB110220'
#f_n     = 'FRB010621'
#f_q     = 'FRB010125'
f_n     = 'FRB110220'
f_q =f_n
#f_q     = 'FRB010621'
f_dir   = '../data/'
#data = np.load(f_dir + 'Crab_52277_022_short.npy')
#data    = np.load(f_dir + 'signal_sim_snr1.npy')
#data    = np.load(f_dir + 'FRB110220.npy')
#data    = np.load(f_dir + 'FRB010125.npy')
data    = np.load(f_dir + f_n + '.npy')
#data    = np.random.random(data.shape)
print data.shape,'**'
#data[50:70,:]=1000
#T = 2
#F = 3000/T
#Value = 5
#data[2400:2400+F,50:50+T] = Value*np.random.random((F,T))

#data[1900:1900+F,170:170+T] = Value*np.random.random((F,T))

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
#data  = data[:,:512]
pixel = 2
area  = (pixel*2)**2
print data.shape,data.max()
freq    = np.load(f_dir + f_q +'_freq.npy')
print freq.shape,'######freq shape#####'
#freq   = np.load(f_dir + 'Crab_52277_022_freq.npy')
#freq    = np.load(f_dir + 'freq.npy')
#freq    = np.load(f_dir + 'FRB110220_freq.npy')
#freq    = np.load(f_dir + 'FRB010125_freq.npy')
#f_n     = 'FRB010125'
fy           = freq**-2
#nbin         = max(data.shape)
nbin        = 512
print data.shape,nbin
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

print data.shape
fft2   = np.fft.fft2(abs(data))
fft2   = np.nan_to_num(fft2)
shift  = np.fft.fftshift(fft2)/np.sqrt(fft2.size)
data   =  shift[1:shift.shape[0]/2+1,shift.shape[1]/2:]
print data.shape
#exit()
#temp =data[-1,0]
data[-5:, :] = 0
data[ :,:5] = 0
#data[-1,0] = temp
#data   = data[:-20,20:]
print data.shape
data2  = data

rang   = data.shape
row    = np.arange(-rang[0]+1,1,dtype=np.float32)
line   = np.arange(rang[1],dtype=np.float32)
angle  = np.nan_to_num(np.arctan(row[:,None]/line[None,:]) / np.pi * 180)
radius = np.sqrt(row[:,None]**2 + line[None,:]**2)
#ang    = angle.reshape(-1)
#rad    = radius.reshape(-1)
ang    = angle.reshape(-1)
rad    = radius.reshape(-1)
ang_min = -5
ang_max = -85
rad_grid = 256#data.shape[0]
ang_grid = 256#data.shape[1]
points = (rad,ang)
data   = data.reshape(-1)
grid_r, grid_a = np.mgrid[0:rad.max():rad_grid*1j, ang_min:ang_max:ang_grid*1j]
polar_matrix   = griddata(points,data,(grid_r,grid_a),method='nearest')
polar_matrix   = np.nan_to_num(polar_matrix)
if rang[0] > rang[1]:
    short =  rang[1]
else:
    short = rang[0]
#polar_matrix[:,-3:] = 0
#polar_matrix[:3,:] = 0
#polar_matrix = polar_matrix[:short,:]
data3 = polar_matrix
#print data3.shape,'***data3.shape'
data4 = abs(data3).sum(axis=0)
#data4[-10:] = data4[:10] = data4[10:-5].mean()

#fit
print data4.mean(),data4.std(),'data4 original mean and std'
tem_data4   = data4
filp=115
prof_data4  = signal.medfilt(data4,filp)
#prof_data4  = data4.mean()
data4      = data4 - prof_data4
tem = data4
data4_max = 0

for i in np.arange(3):
    data4_max = data4_max + data4.max()
    lo  =np.where(data4==data4.max())
    print lo
    data4[lo[0]]=0
data4_std = data4.std()
data4_mean = data4.mean()
data4 = tem
snr  = (data4_max)/data4_std
print 'data4 filt lengt',filp,'SNR',snr, 'mean:',data4_mean, 'old std',data4.std(),'new: std',data4_std
#exit()
#fft1 = np.fft.fft(data3,axis=0)
#fft1 = fft1/np.sqrt(fft1.shape[0])
#shift  = np.fft.fftshift(fft1,axes=0)
#data5 = shift

fft1    = np.fft.fft2(data3)
fft2    = fft2/np.sqrt(fft2.size)
shift   = np.fft.fftshift(fft2)
data5   = shift

data5=abs(data5)
c=data5
print c.max(),c.shape
s = ['rx','r+','ro','r>','r<','rv','r^','rs','r.','r*','bx','b+','bo','b>','b<','bv','b^','sb','b.','b*','gx','g+','go','g>','g<','gv','g^','sg','g.','g*']
temx=np.linspace(ang_min,ang_max,c.shape[1])
temy=np.linspace(-c.shape[0]/2,c.shape[0]/2,c.shape[0])
print temx.shape,temy.shape
d_maxx=0
#plt.figure(figsize=[12,12])
for i in range(area):
    lo    = np.where(data5 == np.max(data5))
  #  print 'i:',i,lo,'Value:',temx[lo[1][0]],temy[lo[0][0]]#data5[lo[0][0],lo[1][0]]
    d_maxx = d_maxx + np.max(data5) 
    data5[lo[0][0],lo[1][0]]=0
#    plt.plot(temx[lo[1][0]],temy[lo[0][0]],s[i],label=str(i))

#plt.legend()
#plt.show()
print c.max(),c.mean(),c.std()
#exit()
data5 = c
lo    = np.where(data5 == np.max(data5))
d_max = 0
data5_n = data5-data5.max()
SNRt = (data5.max()-data5_n.mean())/data5_n.std()
for i in np.arange(-pixel,pixel):
    for j in np.arange(-pixel,pixel):
        d_max += data5[lo[0][0]+i,lo[1][0]+j]
SNR  = d_maxx/((data5[:,-100:-50].mean()+data5[:,-100:-50].std())*area)



plt.figure(figsize=(14,10))
plt.subplot(4,2,1)
plt.title('raw data')
plt.xlabel('time')
plt.ylabel('Frequency')
plt.pcolormesh(abs(data0))
plt.colorbar()

plt.subplot(4,2,2)
plt.title('rebin data')
plt.xlabel('time')
plt.ylabel('wave^2')
plt.pcolormesh(abs(data1))
plt.colorbar()

plt.subplot(4,2,3)
plt.title('2D-FFT ')
plt.xlabel('time axis after FFT')
plt.ylabel('wave^2 axis after FFT')
plt.pcolormesh(abs(data2))
plt.colorbar()

plt.subplot(4,2,4)
plt.title('Polar coordinates transform')
plt.xlabel('angle')
plt.ylabel('radius')
x_axis = np.linspace(ang_min,ang_max,data3.shape[1])
y_axis = np.arange(data3.shape[0])
plt.pcolormesh(x_axis,y_axis,abs(data3))
plt.colorbar()

plt.subplot(4,2,5)
plt.title('sum along radius at each angle(SNR:'+str(int(snr))+')')
plt.xlabel('angle from '+str(ang_min)+' to '+str(ang_max)+' degree')
plt.ylabel('abs value of sum at each angle')
plt.plot(x_axis,data4)
plt.plot(x_axis,tem_data4)
plt.plot(x_axis,prof_data4)

plt.subplot(4,2,6)
data5  = abs(data5)
x_axis = np.linspace(ang_min,ang_max,data5.shape[1])
y_axis = np.arange(-data5.shape[0]/2,data5.shape[0]/2)
plt.xlabel('SNR:'+str(SNR))
plt.pcolormesh(x_axis,y_axis,abs(data5))
plt.colorbar()
lo=np.where(data5==data5.max())
x_max=x_axis[lo[1][0]]
y_max=y_axis[lo[0][0]]
plt.plot(x_max,y_max)#,'ro')


plt.subplot(4,2,7)
data6=abs(data5).sum(axis=0)
#print data6.max(),data6.
#exit()
tem1 = data6
data6_max=0
prof_data6  = signal.medfilt(data6,115)
data6      = data6 - prof_data6
tem = data6
for i in np.arange(2*pixel):
    data6_max = data6_max + data6.max()
    lo=np.where(data6 == data6.max())
    data6[lo[0]]= 0
data6_std  =   data6.std()
data6_mean  = data6.mean()
SNR_s   =   (data6_max)/data6_std
print 'data6',SNR_s, 'new',data6.mean(),'old',tem.mean(),'new std',data6_std,'old std',tem.std()
plt.xlabel('SNR:'+str(SNR_s))
plt.plot(x_axis, tem)
plt.plot(x_axis, prof_data6)
plt.plot(x_axis,tem1)
plt.savefig(f_q)
plt.show()
#print 'F:',F,"and  T:",T, '  Value:',Value
#print 'pixel:',F*T
print 'snr:',snr,'  SNR:',SNR,' SNRt',SNRt,'SNR_s',SNR_s
