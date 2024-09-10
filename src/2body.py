import xymass
import orbit_animation
import numpy as np
import scipy
import astropy.units as u
from astropy.io import fits
import matplotlib.pyplot as plt

period=0.01*u.yr
eccentricity=0.
mass_primary=1.*u.M_sun
mass_secondary=1.*u.M_sun
longitude=(0.*u.deg).to(u.rad) #longitude of ascending node
inclination=(90.*u.deg).to(u.rad) #inclination angle

mass_ratio=(mass_secondary/mass_primary).value
#sample f_period = time / period uniformly over one full period
f_period=np.linspace(0,1,200)

#calculate orbit
sample_orbit=xymass.sample_orbit_2body(f_period,period=period,eccentricity=eccentricity,\
                                       mass_primary=mass_primary,mass_ratio=mass_ratio,\
                                       longitude=longitude,inclination=inclination)

#animate positions of reduced mass, particle 1 and particle 2
#orbit_animation.animation_2body_r(sample_orbit,animation_filename=None)

#animate velocities of reduced mass, particle 1 and particle 2
#orbit_animation.animation_2body_v(sample_orbit,animation_filename=None)


def get_spectrum(wav,flux):
    class spectrum:
        def __init__(self,wav=None,flux=None):
            self.wav=wav
            self.flux=flux
    return spectrum(wav=wav,flux=flux)

wav=fits.open('WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')[0].data
sun=fits.open('lte05800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits')[0].data
#sun=10*wav

npix=10000
gaussian_filter_scale=5.
flux=np.zeros((len(sample_orbit.r_xyz),npix))
flux1=np.zeros((len(sample_orbit.r_xyz),npix))
flux2=np.zeros((len(sample_orbit.r_xyz),npix))
wav_interp=np.linspace(5100,5200,npix)
for i in range(0,len(sample_orbit.r_xyz)):
    wav1=wav*(1.+sample_orbit.v1_obs_xyz.T[2][i].to(u.km/u.s).value/3.0e+5)
    wav2=wav*(1.+sample_orbit.v2_obs_xyz.T[2][i].to(u.km/u.s).value/3.0e+5)
    flux1_0=np.interp(wav_interp,wav1,sun)
    flux2_0=np.interp(wav_interp,wav2,sun)
    flux_0=flux1_0+flux2_0
    flux[i]=scipy.ndimage.gaussian_filter(flux_0,sigma=gaussian_filter_scale)
    flux1[i]=scipy.ndimage.gaussian_filter(flux1_0,sigma=gaussian_filter_scale)
    flux2[i]=scipy.ndimage.gaussian_filter(flux2_0,sigma=gaussian_filter_scale)
    flux[i]/=np.max(flux[i])
    flux1[i]/=np.max(flux1[i])
    flux2[i]/=np.max(flux2[i])
spec1=get_spectrum(wav=wav_interp,flux=flux1)
spec2=get_spectrum(wav=wav_interp,flux=flux2)
spec=get_spectrum(wav=wav_interp,flux=flux)
#spectrum=np.random
orbit_animation.animation_2body_spec(sample_orbit,spec,spec1,spec2,animation_filename=None)
