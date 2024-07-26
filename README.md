# xymass

A package for generating random samples of 2D stellar positions from common surface density models (plummer, exponential, generalized plummer), and random samples of masses from common initial mass functions (Salpeter, Kroupa, broken power law, log-normal), optionally including biary star components, with positions obtained either by some specified separation function or by physical calculation of 2-body orbits.

Author: Matthew G. Walker (2024) 

# Instructions 

* Install xymass. You can either pip install the released version or install from github

```
pip install xymass
```
# Available Models for 2D position

Options for modeling 2D positions are Plummer ('plum'), exponential ('exp'), uniform ('uni') and the projection of an alpha/beta/gamma model with alpha=2 ('a2bg').  

The Plummer model has the form $\Sigma(R)=\frac{\Sigma_0}{(1+R^2/R_p^2)^2}$.

The Exponential model has the form $\Sigma(R)=\Sigma_0\exp[-R/R_e]$.

The uniform model has the form $\Sigma(R)=\Sigma_0$.

The a2bg model has the form $\Sigma(R)=2\int_{R}^{\infty}\frac{r \nu(r)dr}{\sqrt{r^2-R^2}}$, where the 3D profile is $\nu(r)=\frac{\nu_0}{(r/r_s)^{\gamma}[1+r^2/r_s^2]^{(\beta-\gamma)/2}}$.

For all models, the constant $\Sigma_0$ is determined by choice of number of objects and other model parameters (e.g., scale radius and other shape parameters if applicable).

# Available Models for stellar mass function

Options for modeling the stellar IMF are Salpeter ('salpeter'), Kroupa ('kroupa'), broken power law ('bpl') and log-normal ('lognormal').  

The Salpeter model has the form $dN/dM=k M^{-\alpha}$.

The Kroupa model has the form $dN/dM = k_1m^{-\alpha_1}$ for $m\leq m_{\rm break,1}$, $dN/dM=k_2m^{-\alpha_2}$ for $m_1< m\leq m_{\rm break,2}$, $dN/dM=k_3m^{-\alpha_3}$ for $m>m_{\rm break,3}$.

The broken power law model has the form $dN/dM = k_1m^{-\alpha_1}$ for $m\leq m_{\rm break,1}$, $dN/dM=k_2m^{-\alpha_2}$ for $m>m_{\rm break,2}$.

The log-normal model has the form $dN/d\log M = \mathcal{N}(\overline{\log_{10}[M/M_{\odot}]},\sigma_{\log_{10}[M/M_{\odot}})$, where $\mathcal{N}(\overline{x},\sigma_x)$ is the normal distribution with mean $\overline{x}$ and standard deviation $\sigma_x$.

# Usage 

In order to sample 2D positions, specify sample size and analytic model ('plum', 'exp', 'uni', 'a2bg'):

```sample_xy=xymass.sample_r2d(1000,'plum')```

Optionally, specify nonzero ellipticity and position angle (degrees):

```sample_xy=xymass.sample_r2d(1000,'plum',ellipticity=0.4,position_angle=35)```

The model scale radius is 1 by default.  For other values, specify as

```sample_xy=xymass.sample_r2d(1000,'plum',r_scale=1.9,ellipticity=0.4,position_angle=35)```

The returned object contains sampled positions x, y, z (although the 'z' component is zero by construction) 'elliptical' radii (semi-major axis of ellipse centered on origin that passes through sampled position), model parameters, and a function that returns the expected number density at a given elliptical radius.

If using the 'a2bg' model (alpha/beta/gamma model with alpha=2), must specify beta and gamma, e.g.:

```sample_xy=xymass.sample_r2d(1000,'a2bg',r_scale=5.3,beta=5.4,gamma=0.9,ellipticity=0.4,position_angle=35)```


In order to sample stellar masses, specify sample size and analytic model ('salpeter', 'kroupa', 'lognormal', 'bpl'):

 ```sample_mass=xymass.sample_imf(1000,'kroupa')```

The returned object contains sampled masses, model parameters (including normalization constants), and a function that returns the expected probability density, dN/dmass (normalized so that integral over dN/dmass between m_min and m_max is 1), at a given mass value.  

Default values of model parameters are given below.  Different values can be specified when calling, e.g.

 ```sample_mass=xymass.sample_imf(1000,'kroupa',alpha1=0.5,alpha2=1.5,alpha3=2.5,m1_break=0.1,m2_break=0.6)```

Default values:

Salpeter: alpha=2.3

Kroupa: alpha1=0.3, alpha2=1.3, alpha3=2.3, m1_break=0.08, m2_break=0.5

broken power law: alpha1=1.3, alpha2=2.3, m_break=0.5

log-normal: mean=0.08, std=0.7 

In order to replace some fraction f_binary of sampled coordinates with positions of individual members of binary star systems based on 2body orbit calculations:

 ```sample_xy_withbinaries=xymass.add_binaries_physical(sample_xy.r_xyz,sample_mass.mass*u.M_sun,f_binary=f_binary,m_min=m_min,binary_model='Raghavan2010')```

The positions of the stars within the binary systems are sampled by drawing from a 2-body Keplerian orbit (randomly-drawn physical and observational parameters), with binary orbital parameters specified either by the user (see examples) or,as here, according to the parameters inferred by Raghavan et al. 2010 (Duquennoy & Mayor 1991 also available as 'DM91).  m_min should be set to the minimum mass allowed for the secondary (e.g., hydrogen-burning limit).  

In order to sample phase-space coordinates (with respect to center of mass) of binary companions directly:

 ```sample_orbit=xymass.sample_orbit_2body(f_period,period=period,eccentricity=eccentricity,mass_primary=mass_primary,mass_ratio=mass_ratio,longitude=longitude,inclination=inclination)```

where f_period specifies the time of observation as a fraction of the orbital period, and the other parameters specify the orbital and observational parameters (all except f_period, eccentricity and mass ratio must have units specified via astropy.units).  All can be scalars or numpy arrays.   
 
# Examples 

For examples of sampling 2D positions, see the [notebook](examples/sample_r2d.ipynb) in the examples folder.

For examples of sampling initial masses, see the [notebook](examples/sample_imf.ipynb) in the examples folder.

For examples of sampling position and velocity components due to binary star orbits, see the [notebook](examples/sample_orbit.ipynb) in the examples folder.

For examples of adding physically-calculated binary star positions to the sampled 2D positions, see the [notebook](examples/sample_r2d_with_binaries.ipynb) in the examples folder.

For examples of generating a synthetic stellar population, optionally with binary stars, with specific filter magnitudes calculated using the (MIST) isochrones package, including blending effects due to instrumental angular resolution, see the [notebook](examples/sample_r2d_with_binaries_blend.ipynb) in the examples folder.

# Acknowledgement

