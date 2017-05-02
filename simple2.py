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

#f_n      = 'signal_sim_snr_1_c'
f_n	= 'filtered_short'

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
data	= data[:,0,:]
#print data.shape,'2222'
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
freq    = np.load(f_dir +'freq.npy')
#freq   = np.load(f_dir + 'Crab_52277_022_freq.npy')
#freq    = np.load(f_dir + 'freq.npy')
#freq    = np.load(f_dir + 'FRB110220_freq.npy')
#freq    = np.load(f_dir + 'FRB010125_freq.npy')
#f_n     = 'FRB010125'
fy           = freq**-2
#nbin         = max(data.shape)
#nbin        = data.shape[1]
#print data.shape,nbin
	#nbin         = 1024
nbin = np.arange(data.shape[0]/4,data.shape[0]*6)
#nbin = np.arange(1900,2100)
#nbin = np.arange(1021,1024)
print 'Length of nbin:',len(nbin)
a=[]
e=[]
b=[]
d=[]
f=[]
for i in np.arange(len(nbin)):
	print 'i step: ',i
	b.append(nbin[i])
	bins	= nbin[i]
	bw	= 200./bins
	f_axis       = np.linspace(fy.min(),fy.max(),bins)
	b_aray,b_edg = np.histogram(fy,bins=bins)
	
	data0 = data 
	data_0nan = np.nan_to_num(data)
	data1 = np.zeros((bins,data.shape[1]),data.dtype)
	for ii in np.arange(bins):
	    length  =  b_aray[ii]
	    if  ii == 0:
	        index = 0
	    elif ii == 1:
	        index = b_aray[0]
	    else:
	        index = b_aray[0:ii].sum()
	    for iii in np.arange(length):
	        data1[ii,:]+=data_0nan[iii+index,:]
	tem = np.ones(data.shape[1])
	bin_weight = b_aray[:,None] * tem[None,:]
	data1 = data1 / 1.0 / bin_weight
	data1 = np.nan_to_num(data1)
	data  = data1
	
#	print data.shape
	fft2   = np.fft.fft2(data)
	fft2   = np.nan_to_num(fft2)
	shift  = np.fft.fftshift(fft2)/np.sqrt(fft2.size)
	data   =  shift[1:shift.shape[0]/2+1,shift.shape[1]/2:]
	#data   =  shift[:shift.shape[0]/2,shift.shape[1]/2:]
	#exit()
	data[-2:, :] = 0
	data[ :,:2] = 0
	#data   = data[:,1:]
	
	
#	print data.shape,'/*/*/*/*/*'
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
	rad_grid = data.shape[0]*1
	ang_grid = data.shape[0]*1
	points = (rad,ang)
	data   = data.reshape(-1)
	grid_r, grid_a = np.mgrid[0:rang[1]:rad_grid*1j, ang_min:ang_max:ang_grid*1j]
	polar_matrix   = griddata(points,data,(grid_r,grid_a),method='nearest')
	polar_matrix   = np.nan_to_num(polar_matrix)
	#print rang,rad.max(),'1231',polar_matrix.shape
	if rang[0] > rang[1]:
	    short =  rang[1]
	else:
	    short = rang[0]
	#polar_matrix[:,-3:] = 0
	#polar_matrix[:3,:] = 0
	polar_matrix = polar_matrix[:short,:]
	data3 = polar_matrix
	#print data3.shape,'***data3.shape'
	data4 = abs(data3).sum(axis=0)
	#data4[-10:] = data4[:10] = data4[10:-5].mean()
	
	#fit
	prof_data  = signal.medfilt(data4,15)
	data4      = data4 - prof_data
	
	
#	fft1 = np.fft.fft(data3,axis=0)
#	fft1 = fft1/np.sqrt(fft1.shape[0])
##	shift  = np.fft.fftshift(fft1,axes=0)
#	data5 = shift
	
	
	snr  	= (data4.max()-data4.mean())/data4[100:].std()
	snr_o	= (data4.max()-data4.mean())/data4.std()
#	data5=abs(data5)
#	c=data5
#	for i in range(15):
#	    lo    = np.where(data5 == np.max(data5))
	#    print 'i:',i,lo,'Value:',data5[lo[0][0],lo[1][0]]
#	    data5[lo[0][0],lo[1][0]]=0
	#    plt.plot(lo[0][0],lo[1][0],'o',label=str(i))
	#    plt.legend()
	#plt.show()
#	print c.max(),c.mean(),c.std(),'asdfasdf'
	#exit()
#	data5 = c
#	lo    = np.where(data5 == np.max(data5))
#	d_max = 0
#	SNRt = (data5.max()-data5.mean())/data5.std()
#	for i in np.arange(-pixel,pixel):
#	    for j in np.arange(-pixel,pixel):
#	        d_max += data5[lo[0][0]+i,lo[1][0]+j]
#	SNR  = (d_max)/((data5[:,-100:].std()+data5[:,:100].mean())*area)
	
	data=data0
#	print data4.max() , snr,'*****', (data4[50:]).std(),'\n\n'
#	print a[i],b[i],e[i],'a,b,e\n\n'
	a.append(snr)
	e.append(bw)
	f.append(snr_o)
	d.append(data4.max())
	print 'nbin: ',bins, '  SNR:',snr,'------ SNR_O: ',snr_o
	
	print 'defference between snr and snr_o',snr-snr_o
	print data4.shape,'%%%%% std defference of std([:]-[100:)]:',data4.std()-data4[100:].std()
	print '****\n\n'
snr	=	np.array(a)
nbin	=	np.array(b)
maxv	=	np.array(d)
bw	=	np.array(e)
snr_o	=	np.array(f)
np.save('snr',snr)
np.save('nbin',nbin)
np.save('maxvalue',maxv)
np.save('bw',bw)
np.save('snr_o',snr_o)
#plt.plot(b,a)
#plt.grid()
#plt.show()

intsnr	= signal.medfilt(snr,55)
intmax  = signal.medfilt(maxv,55)

plt.figure(figsize=(14,14))

plt.subplot(2,2,1)
plt.title('snr-bw')
plt.plot(bw,snr,label='SNR')
plt.plot(bw,maxv,label='maxv')
plt.plot(bw,intsnr,label='medfit of SNR')
plt.xlabel('Channel width after rebin'); plt.legend(loc='best')

plt.subplot(2,2,2)
plt.title('max value - bw')
plt.plot(bw,maxv,label='maxv')
plt.plot(bw,intmax,label='medfit of max value')
plt.xlabel('Channel width after rebin'); plt.legend(loc='best')

plt.subplot(2,2,3)
plt.title('snr-nbin')
plt.plot(nbin,snr,label='SNR')
plt.plot(nbin,maxv,label='maxv')
plt.plot(nbin,intsnr,label='medfit of SNR')
plt.xlabel('Number of bins after rebin'); plt.legend(loc='best')

plt.subplot(2,2,4)
plt.title('max value-nbin')
plt.plot(nbin,maxv,label='maxv')
plt.plot(nbin,intmax,label='medfit of max value')
plt.xlabel('Number of bins after rebin'); plt.legend(loc='best')

plt.show()
exit()
	
#	plt.figure(figsize=(21,14))
#	plt.subplot(3,2,1)
#	plt.title('raw data')
#	plt.xlabel('time')
#	plt.ylabel('Frequency')
#	plt.pcolormesh(abs(data0))
#	plt.colorbar()
#	
#	plt.subplot(3,2,2)
#	plt.title('rebin data')
#	plt.xlabel('time')
#	plt.ylabel('wave^2')
#	plt.pcolormesh(abs(data1))
#	plt.colorbar()
#	
#	plt.subplot(3,2,3)
#	plt.title('2D-FFT ')
#	plt.xlabel('time axis after FFT')
#	plt.ylabel('wave^2 axis after FFT')
#	plt.pcolormesh(abs(data2))
#	plt.colorbar()
#	
#	plt.subplot(3,2,4)
#	plt.title('Polar coordinates transform')
#	plt.xlabel('angle')
#	plt.ylabel('radius')
#	x_axis = np.linspace(ang_min,ang_max,data3.shape[1])
#	y_axis = np.arange(data3.shape[0])
#	plt.pcolormesh(x_axis,y_axis,abs(data3))
#	plt.colorbar()
#	
#	plt.subplot(3,2,5)
#	plt.title('sum along radius at each angle(SNR:'+str(int(snr))+')')
#	plt.xlabel('angle from '+str(ang_min)+' to '+str(ang_max)+' degree')
#	plt.ylabel('abs value of sum at each angle')
#	plt.plot(x_axis,data4,'b',label='1st FFT')
#	p=abs(data5).sum(axis=0)
#	plt.plot(x_axis,p,'r',label='2nd FFT')
#	plt.legend(loc='best')
#	
#	plt.subplot(3,2,6)
#	x_axis = np.linspace(ang_min,ang_max,data5.shape[1])
#	y_axis = np.arange(-data5.shape[0]/2,data5.shape[0]/2)
#	plt.xlabel('SNR:'+str(SNR))
#	#p=abs(data5).sum(axis=0)
#	#plt.plot(x_axis,p)
#	plt.pcolormesh(x_axis,y_axis,abs(data5))
#	plt.colorbar()
#	plt.savefig(f_n)
#	plt.show()
#	#print 'F:',F,"and  T:",T, '  Value:',Value
#	#print 'pixel:',F*T
#	print 'snr:',snr,'  SNR:',SNR,' SNRt',SNRt
