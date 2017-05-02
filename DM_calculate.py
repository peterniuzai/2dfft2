import numpy as np
import sys
def DM_calculate(f_axis,t_rsl,location,nbin):
    d     = 4.148908e6 # (ms)
    d     = d / 1000
    
    loc   = np.float32(90 - location)
    ang   = np.tan(loc/180*np.pi)
    f_rsl = (f_axis[-1] - f_axis[0])/nbin
    unit  = f_rsl / t_rsl
    k     = ang * unit
    DM    = 1/k/d
    return DM

def degree_calculate(f_axis,t_rsl,DM,nbin):
    d      = 4.148908e6 # (ms)
    d      = d / 1000

    f_rsl  = (f_axis[-1] - f_axis[0])/nbin
    unit   = f_rsl / t_rsl
    k      = 1 / DM / d
    ang    = k / unit
    degree = np.degrees( np.arctan(ang) )     
    degree = 90 - degree
    return degree

if __name__ == '__main__':
    ang     = np.float(sys.argv[1])
    DM      = np.float(sys.argv[2])
    freq    = np.load('/home/ycli/data/burst_data/freq.npy')
    t_rsl   = 0.001024    #second
    f_rsl   = 0.048828125 #Mhz
    f_axis  = freq ** -2
    nbin    = 1024
    dm = DM_calculate(f_axis,t_rsl,ang,nbin)
    degree = degree_calculate(f_axis,t_rsl,DM,nbin)
    print 'the DM is : ',dm ,'pc*cm^-3'
    print 'The degree is :',degree, ' degree'
