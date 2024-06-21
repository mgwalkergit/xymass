# xymass

A package for generating random samples of 2D stellar positions from common surface density models (plummer, exponential, generalized plummer), and random samples of masses from common initial mass functions (Salpeter, Kroupa, broken power law, log-normal).

Author: Matthew G. Walker (2024) 

# Instructions 

* Install xymass. You can either pip install the released version or install from github

```
pip install xymass
```
# Available Models for 2D position

Options for modeling 2D positions are Plummer ('plum'), exponential ('exp'), uniform ('uni') and the projection of an alpha/beta/gamma model with alpha=2 ('2bg').  

The Plummer model has the form $\Sigma(R)=\frac{\Sigma_0}{(1+R^2/R_p^2)^2}$.

The Exponential model has the form $\Sigma(R)=\Sigma_0\exp[-R/R_e]$.

The uniform model has the form $\Sigma(R)=\Sigma_0$.

The 2bg model has the form $\Sigma(R)=2\int_{R}^{\infty}\frac{r \nu(r)dr}{\sqrt{r^2-R^2}}$, where the 3D profile is $\nu(r)=\frac{\nu_0}{(r/r_s)^{\gamma}[1+r^2/r_s^2]^{(\beta-\gamma)/2}}$.

# Available Models for initial mass function

Options for modeling the IMF are Salpeter ('salpeter'), Kroupa ('kroupa'), broken power law ('bpl') and log-normal ('lognormal').  

The Salpeter model has the form $dN/dM=k M^{-\alpha}$.

The Kroupa model has the form $dN/dM = k_1m^{-\alpha_1}$ for $m\leq m_{\rm break,1}$, $dN/dM=k_2m^{-\alpha_2}$ for $m_1< m\leq m_{\rm break,2}$, $dN/dM=k_3m^{-\alpha_3}$ for $m>m_{\rm break,3}$.

The broken power law model has the form $dN/dM = k_1m^{-\alpha_1}$ for $m\leq m_{\rm break,1}$, $dN/dM=k_2m^{-\alpha_2}$ for $m>m_{\rm break,2}$.

The log-normal model has the form $dN/d\log M = \mathcal{N}(\overline{\log_{10}[M/M_{\odot}]},\sigma_{\log_{10}[M/M_{\odot}})$, where $\mathcal{N}(\overline{x},\sigma_x)$ is the normal distribution with mean $\overline{x}$ and standard deviation $\sigma_x$.

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

Salpeter: alpha=2.3

Kroupa: alpha1=0.3, alpha2=1.3, alpha3=2.3, m1_break=0.08, m2_break=0.5

broken power law: alpha1=1.3, alpha2=2.3, m_break=0.5

log-normal: mean=0.08, std=0.7 

# Examples 

For examples of sampling 2D positions, see the [notebook](examples/xymass_sample_r2d_example.ipynb) in the examples folder.

For examples of sampling initial masses, see the [notebook](examples/xymass_sample_imf_example.ipynb) in the examples folder.

# Acknowledgement
