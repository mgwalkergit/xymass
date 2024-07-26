import xymass
import numpy as np
import scipy
import astropy.units as u
import matplotlib.pyplot as plt

n_object=100000 #sample size
f_binary=0.5 #binary fraction

m_min=0.1
mass_primary=xymass.sample_imf(size=n_object,model='kroupa',m_min=m_min).mass*u.M_sun
mass_secondary=xymass.sample_imf(size=n_object,model='kroupa',m_min=m_min).mass*u.M_sun

r2d=xymass.sample_r2d(size=n_object,model='plum',r_scale=100.*u.pc)
r2d_with_binaries_raghavan=xymass.add_binaries_physical(r2d.r_xyz,mass_primary,f_binary=f_binary,m_min=m_min,binary_model='Raghavan2010')

mass_ratio=np.random.uniform(size=n_object,low=m_min/mass_primary.value,high=1.) #mass ratio M_secondary / M_primary
s_min=0.1*u.AU
s_max=10000*u.AU
bpl_s_break=1000.*u.AU
bpl_alpha1=1.01
bpl_alpha2=4.
pl_alpha=3.001
lognormal_loc=100.*u.AU
lognormal_scale=10.*u.AU

r2d_with_binaries_bpl=xymass.add_binaries_func(r2d.r_xyz,separation_func='bpl',mass_primary=mass_primary,f_binary=f_binary,mass_ratio=mass_ratio,s_min=s_min,s_max=s_max,s_break=bpl_s_break,alpha1=bpl_alpha1,alpha2=bpl_alpha2,projected=True)
r2d_with_binaries_lognormal=xymass.add_binaries_func(r2d.r_xyz,separation_func='lognormal',mass_primary=mass_primary,mass_secondary=mass_secondary,f_binary=f_binary,s_min=s_min,s_max=s_max,loc=lognormal_loc,scale=lognormal_scale,projected=True)
r2d_with_binaries_pl=xymass.add_binaries_func(r2d.r_xyz,separation_func='pl',mass_primary=mass_primary,mass_secondary=mass_secondary,f_binary=f_binary,s_min=s_min,s_max=s_max,alpha=pl_alpha,projected=True)
