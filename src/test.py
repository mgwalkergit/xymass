import xymass
import numpy as np
import scipy
import astropy.units as u
import matplotlib.pyplot as plt

n_object=10000 #sample size
f_binary=0.5 #binary fraction

m_min=0.1
mass_primary=xymass.sample_imf(size=n_object,model='kroupa',m_min=m_min).mass*u.M_sun
mass_secondary=xymass.sample_imf(size=n_object,model='kroupa',m_min=m_min).mass*u.M_sun 

r2d=xymass.sample_r2d(size=n_object,model='plum',r_scale=100.*u.pc,ellipticity=0.,position_angle=0.)

s_min=1.*u.AU #minimum separation
s_max=100000*u.AU #maximum separation

r2d_with_binaries_opik=xymass.add_binaries_func(r2d.r_xyz,\
                                              separation_func='opik',mass_primary=mass_primary,\
                                              mass_secondary=mass_secondary,f_binary=f_binary,s_min=s_min,\
                                              s_max=s_max,projected=True) #project=True if separation function is 2D, False if 3D

