# xymass

A package for generating random samples of 2D stellar positions from common surface density models (plummer, exponential, generalized plummer), and random samples of masses from common initial mass functions (Salpeter, Kroupa, broken power law, log-normal).

Author: Matthew G. Walker (2024) 

# Instructions 

* Install xymass. You can either pip install the released version or install from github

```
pip install xymass
```
# Usage 

In order to sample 2D positions, specify sample size and analytic model ('plum', 'exp', 'uni', '2bg'):

```sample_xy=xymass.sample_r2d(1000,'plum')```

Optionally, specify nonzero ellipticity and position angle (degrees):

```sample_xy=xymass.sample_r2d(1000,'plum',ellipticity=0.4,position_angle=35)```

The model scale radius is 1 by default.  For other values, specify as

```sample_xy=xymass.sample_r2d(1000,'plum',r_scale=1.9,ellipticity=0.4,position_angle=35)```

The returned object contains sampled positions x, y, 'elliptical' radii (semi-major axis of ellipse centered on origin that passes through sampled position), model parameters, and a function that returns the expected number density at a given elliptical radius.

If using the '2bg' model (alpha/beta/gamma model with alpha=2), must specify beta and gamma, e.g.:

```sample_xy=xymass.sample_r2d(1000,'2bg',r_scale=5.3,beta=5.4,gamma=0.9,ellipticity=0.4,position_angle=35)```


 In order to sample stellar initial masses, specify sample size and analytic model ('salpeter', 'kroupa', 'lognormal', 'bpl'):

 ```sample_mass=xymass.sample_imf(1000,'kroupa')```

The returned object contains sampled masses, model parameters (including normalization constants), and a function that returns the expected probability density, dN/dmass (normalized so that integral over dN/dmass between m_min and m_max is 1), at a given mass value.

Default values of model parameters are given below.  Different values can be specified when calling, e.g.

 ```sample_mass=xymass.sample_imf(1000,'kroupa',alpha1=0.5,alpha2=1.5,alpha3=2.5,m1_break=0.1,m2_break=0.6)```

Default values:

salpeter_alpha=2.3

kroupa_alpha1=0.3

kroupa_alpha2=1.3

kroupa_alpha3=2.3

kroupa_m1_break=0.08 (Msun)

kroupa_m2_break=0.5 (Msun)

bpl_alpha1=1.3

bpl_alpha2=2.3

bpl_m_break=0.5 (Msun)

lognormal_mean=0.08 (units of Msun) 

lognormal_std=0.7 (units of Msun)



# Examples 

For examples of sampling 2D positions, see the [notebook](examples/xymass_sample_r2d_example.ipynb) in the examples folder.

For examples of sampling initial masses, see the [notebook](examples/xymass_sample_imf_example.ipynb) in the examples folder.

# Acknowledgement
