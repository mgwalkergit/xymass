import numpy as np
import scipy
import scipy.optimize
import scipy.special
import warnings
from scipy.spatial.transform import Rotation
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from astropy import units as u
import time

def sample_r2d(size,model,**params):#samples from flattened plummer, exponential, or (not flattened) uniform 2d distributions
    
    class r2d:
        def __init__(self,r_ell=None,x=None,y=None,r_xyz=None,ellipticity=None,position_angle=None,r_scale=None,model=None,alpha=None,beta=None,gamma=None,func=None,rhalf_2d=None):
            self.r_ell=r_ell
            self.x=x
            self.y=y
            self.r_xyz=r_xyz
            self.ellipticity=ellipticity
            self.position_angle=position_angle
            self.r_scale=r_scale
            self.model=model
            self.alpha=alpha
            self.beta=beta
            self.gamma=gamma
            self.func=func
            self.rhalf_2d=rhalf_2d

    def flatten_2d(size,params):#computes x,y coordinates (units of r_scale) given ellipticity and position angle (units of R/r_scale**2)
        phi=2.*np.pi*np.random.uniform(low=0.,high=1.,size=size)#azimuthal angle in circular coordinates
        x0,y0=np.cos(phi)*(1.-params['ellipticity']),np.sin(phi)#stretch along x axis
        xflat=x0*np.cos(-params['position_angle']*np.pi/180.)-y0*np.sin(-params['position_angle']*np.pi/180.)#now rotate axes by position angle
        yflat=y0*np.cos(-params['position_angle']*np.pi/180.)+x0*np.sin(-params['position_angle']*np.pi/180.)
        return xflat,yflat

    if not 'r_scale' in params:
        params['r_scale']=1.
        warnings.warn('r_scale not specified, assuming r_scale=1')
    if params['r_scale']<0:
        raise ValueError('r_scale = '+str(params['r_scale'])+' is invalid value, must have r_scale >=0.')
    if (('position_angle' in params)&(not 'ellipticity' in params)):
        raise ValueError('specified position_angle but not ellipticity')
    if (('position_angle' not in params)&('ellipticity' in params)):
        raise ValueError('specified ellipticity but not position_angle')        
    if ((not 'ellipticity' in params)&(not 'position_angle' in params)):
        params['ellipticity']=0.
        params['position_angle']=0.
        if not model=='uni':
            warnings.warn('ellipticity and position_angle not specified, assuming ellipticity=0')
    if ((params['ellipticity']<0.)|(params['ellipticity']>1.)):
        raise ValueError('ellipticity = '+str(params['ellipticity'])+' is invalid value, must be between 0 and 1')
    if ((model=='uni')&(params['ellipticity']!=0)):
        warnings.warn('specified uniform distribution with nonzero ellipticity!')
    if model=='2bg':
        if 'beta' not in params:
            raise ValueError('must specify beta and gamma for 2bg model')
        if 'gamma' not in params:
            raise ValueError('must specify beta and gamma for 2bg model')
        
    flat_x,flat_y=flatten_2d(size,params)
    uni=np.random.uniform(low=0.,high=1.,size=size)
    
    if model=='plum':
        
        bigsigma0=size/np.pi/params['r_scale']**2
        rhalf_2d=params['r_scale']
        
        def func(x):
            return bigsigma0/(1+x**2)**2
        r=np.sqrt(uni/(1.-uni))#elliptical radius
        return r2d(r_ell=r*params['r_scale'],x=r*flat_x*params['r_scale'],y=r*flat_y*params['r_scale'],r_xyz=np.c_[r*flat_x*params['r_scale'],r*flat_y*params['r_scale'],np.zeros(len(r),dtype=float)],ellipticity=params['ellipticity'],position_angle=params['position_angle'],r_scale=params['r_scale'],model=model,func=func,rhalf_2d=rhalf_2d)

    if model=='exp':
        
        bigsigma0=size/2/np.pi/params['r_scale']**2
        rhalf_2d=1.67835*params['r_scale']

        def func(x):
            return bigsigma0*np.exp(-x)
        def findbigr_exp(x,uni):
            return 1.-(1.+x)*np.exp(-x)-uni
        low0=0.
        high0=1.e10
        r=[]
        for i in range(0,len(uni)):
            r.append(scipy.optimize.brentq(findbigr_exp,low0,high0,args=uni[i],xtol=1.e-12,rtol=1.e-6,maxiter=100,full_output=False,disp=True))#elliptical radius
        r=np.array(r)
        return r2d(r_ell=r*params['r_scale'],x=r*flat_x*params['r_scale'],y=r*flat_y*params['r_scale'],r_xyz=np.c_[r*flat_x*params['r_scale'],r*flat_y*params['r_scale'],np.zeros(len(r),dtype=float)],ellipticity=params['ellipticity'],position_angle=params['position_angle'],r_scale=params['r_scale'],model=model,func=func,rhalf_2d=rhalf_2d)

    if model=='2bg':
        
        bigsigma0=size*(params['beta']-3)*scipy.special.gamma((params['beta']-params['gamma'])/2)/4/np.sqrt(np.pi)/scipy.special.gamma((3-params['gamma'])/2)/scipy.special.gamma(params['beta']/2)/params['r_scale']**2

        def rootfind_2bg_2d(x,beta,gamma):
            return 0.5-np.sqrt(np.pi)*scipy.special.gamma((beta-gamma)/2)/2/scipy.special.gamma(beta/2)/scipy.special.gamma((3-gamma)/2)*x**(3-beta)*scipy.special.hyp2f1((beta-3)/2,(beta-gamma)/2,beta/2,-1/x**2)

        low0=1.e-10
        high0=1.e+10
        rhalf_2d=params['r_scale']*scipy.optimize.brentq(rootfind_2bg_2d,low0,high0,args=(params['beta'],params['gamma']),xtol=1.e-12,rtol=1.e-6,maxiter=1000,full_output=False,disp=True)

        def func(x):
            return bigsigma0*x**(1-params['beta'])*scipy.special.hyp2f1((params['beta']-1)/2,(params['beta']-params['gamma'])/2,params['beta']/2,-1/x**2)            
        def findbigr_2bg(x,uni,beta,gamma):
            return 1-np.sqrt(np.pi)/2*scipy.special.gamma((beta-gamma)/2)/scipy.special.gamma(beta/2)/scipy.special.gamma((3-gamma)/2)*x**(3-beta)*scipy.special.hyp2f1((beta-3)/2,(beta-gamma)/2,beta/2,-1/x**2)-uni
        low0=1.e-30
        high0=1.e10
        r=[]
        for i in range(0,len(uni)):
            r.append(scipy.optimize.brentq(findbigr_2bg,low0,high0,args=(uni[i],params['beta'],params['gamma']),xtol=1.e-12,rtol=1.e-6,maxiter=100,full_output=False,disp=True))#eliptical radius
        r=np.array(r)
        return r2d(r_ell=r*params['r_scale'],x=r*flat_x*params['r_scale'],y=r*flat_y*params['r_scale'],r_xyz=np.c_[r*flat_x*params['r_scale'],r*flat_y*params['r_scale'],np.zeros(len(r),dtype=float)],ellipticity=params['ellipticity'],position_angle=params['position_angle'],r_scale=params['r_scale'],model=model,beta=params['beta'],gamma=params['gamma'],func=func,rhalf_2d=rhalf_2d)
    
    if model=='uni':
        bigsigma0=size/np.pi/params['r_scale']**2
        rhalf_2d=np.sqrt(0.5)*params['r_scale']
        
        def func(x):
            return bigsigma0*x/x
        r=np.sqrt(uni)#elliptical radius (can in practice be elliptical if nonzero ellipticity is specified)
        return r2d(r_ell=r*params['r_scale'],x=r*flat_x*params['r_scale'],y=r*flat_y*params['r_scale'],r_xyz=np.c_[r*flat_x*params['r_scale'],r*flat_y*params['r_scale'],np.zeros(len(r),dtype=float)],ellipticity=params['ellipticity'],position_angle=params['position_angle'],r_scale=params['r_scale'],model=model,func=func,rhalf_2d=rhalf_2d)

def sample_imf(size,model,**params):
    class imf:
        def __init__(self,model=None,mass=None,mean=None,std=None,alpha=None,alpha1=None,alpha2=None,alpha3=None,m_break=None,m1_break=None,m2_break=None,m_min=None,m_max=None,k=None,k1=None,k2=None,k3=None,func=None):
            self.model=model
            self.mass=mass
            self.mean=mean
            self.std=std
            self.alpha=alpha
            self.alpha1=alpha1
            self.alpha2=alpha2
            self.alpha3=alpha3
            self.m_break=m_break
            self.m1_break=m1_break
            self.m2_break=m2_break
            self.m_min=m_min
            self.m_max=m_max
            self.k=k
            self.k1=k1
            self.k2=k2
            self.k3=k3
            self.func=func

    if not 'm_min' in params:
        params['m_min']=0.1
    if not 'm_max' in params:
        params['m_max']=150.
    if params['m_min']>params['m_max']:
        raise ValueError ('m_min cannot be larger than m_max')
        
    ran1=np.random.uniform(low=0.,high=1.,size=size)
    
    if model=='salpeter':

        if not 'alpha' in params:
            params['alpha']=2.3
        
        k_salpeter=(1.-params['alpha'])/(params['m_max']**(1.-params['alpha'])-params['m_min']**(1.-params['alpha']))
        def salpeter_func(x):
            return k_salpeter*x**-params['alpha']
            
        mass=(params['m_min']**(1.-params['alpha'])+ran1*(params['m_max']**(1.-params['alpha'])-params['m_min']**(1.-params['alpha'])))**(1./(1.-params['alpha']))
        return imf(model=model,mass=mass,alpha=params['alpha'],k=k_salpeter,m_min=params['m_min'],m_max=params['m_max'],func=salpeter_func)

    if model=='lognormal':

        if not 'mean' in params:
            params['mean']=0.08
        if not 'std' in params:
            params['std']=0.7
            
        erf1=scipy.special.erf((np.log10(params['mean'])*np.log(10.)-np.log(params['m_min']))/np.sqrt(2.)/np.log(10.)/params['std'])
        erf2=scipy.special.erf((np.log10(params['mean'])*np.log(10.)-np.log(params['m_max']))/np.sqrt(2.)/np.log(10.)/params['std'])
        k_lognormal=np.sqrt(2./np.pi)/params['std']/(erf1-erf2)
        
        def lognormal_func(x):
            return k_lognormal/x/np.log(10.)*np.exp(-(np.log10(x)-np.log10(params['mean']))**2/2./params['std']**2)
            
        ntotnorm=scipy.special.erf((np.log10(params['mean'])*np.log(10.)-np.log(params['m_min']))/np.sqrt(2.)/np.log(10.)/params['std'])-scipy.special.erf((np.log10(params['mean'])*np.log(10.)-np.log(params['m_max']))/np.sqrt(2.)/np.log(10.)/params['std'])
        erf=scipy.special.erf((np.log10(params['mean'])*np.log(10.)-np.log(params['m_min']))/np.sqrt(2.)/np.log(10.)/params['std'])-ran1*ntotnorm
        mass=np.exp(np.log10(params['mean'])*np.log(10.)-np.sqrt(2.)*np.log(10.)*params['std']*scipy.special.erfinv(erf))
        return imf(model=model,mass=mass,mean=params['mean'],std=params['std'],k=k_lognormal,m_min=params['m_min'],m_max=params['m_max'],func=lognormal_func)
        
    if model=='kroupa':#sample from kroupa IMF, 3 separate power laws with indices -alpha1, -alpha2, -alpha3, break masses at m1_break and m2_break

        if not 'alpha1' in params:
            params['alpha1']=0.3
        if not 'alpha2' in params:
            params['alpha2']=1.3
        if not 'alpha3' in params:
            params['alpha3']=2.3
        if not 'm1_break' in params:
            params['m1_break']=0.08
        if not 'm2_break' in params:
            params['m2_break']=0.5
            
        if params['m1_break']>params['m2_break']:
            raise ValueError ('Kroupa IMF: m1_break cannot be larger than m2_break')
        
        mass=[]

        #get normalization constant for each of three pieces
        k2_over_k1=params['m1_break']**(params['alpha2']-params['alpha1'])
        k3_over_k2=params['m2_break']**(params['alpha3']-params['alpha2'])
        
        if params['m_min']<params['m1_break']:
            
            if params['m_max']>params['m2_break']:
                
                piece1=(params['m1_break']**(1.-params['alpha1'])-params['m_min']**(1.-params['alpha1']))/(1.-params['alpha1'])
                piece2=(params['m2_break']**(1.-params['alpha2'])-params['m1_break']**(1.-params['alpha2']))/(1.-params['alpha2'])
                piece3=(params['m_max']**(1.-params['alpha3'])-params['m2_break']**(1.-params['alpha3']))/(1.-params['alpha3'])

                m0_2=params['m1_break']
                m0_3=params['m2_break']
                                
            if ((params['m1_break']<=params['m_max'])&(params['m_max']<=params['m2_break'])):
                
                piece1=(params['m1_break']**(1.-params['alpha1'])-params['m_min']**(1.-params['alpha1']))/(1.-params['alpha1'])
                piece2=(params['m_max']**(1.-params['alpha2'])-params['m1_break']**(1.-params['alpha2']))/(1.-params['alpha2'])
                piece3=0.

                m0_2=params['m1_break']
                m0_3=params['m2_break']
                                
            if params['m1_break']>params['m_max']:
                
                piece1=(params['m_max']**(1.-params['alpha1'])-params['m_min']**(1.-params['alpha1']))/(1.-params['alpha1'])
                piece2=0.
                piece3=0.
                
                m0_2=params['m1_break']
                m0_3=params['m2_break']
                
        if ((params['m1_break']<=params['m_min'])&(params['m_min']<=params['m2_break'])):
            
            if params['m_max']>params['m2_break']:
                
                piece1=0.
                piece2=(params['m2_break']**(1.-params['alpha2'])-params['m_min']**(1.-params['alpha2']))/(1.-params['alpha2'])
                piece3=(params['m_max']**(1.-params['alpha3'])-params['m2_break']**(1.-params['alpha3']))/(1.-params['alpha3'])
                
                m0_2=params['m_min']
                m0_3=params['m2_break']
                
            if ((params['m1_break']<=params['m_max'])&(params['m_max']<=params['m2_break'])):
                
                piece1=0.
                piece2=(params['m_max']**(1.-params['alpha2'])-params['m_min']**(1.-params['alpha2']))/(1.-params['alpha2'])
                piece3=0.
                
                m0_2=params['m_min']
                m0_3=params['m2_break']
                
        if params['m_min']>params['m2_break']:
            
            if params['m_max']>params['m2_break']:
                
                piece1=0.
                piece2=0.
                piece3=(params['m_max']**(1.-params['alpha3'])-params['m_min']**(1.-params['alpha3']))/(1.-params['alpha3'])
                
                m0_2=params['m1_break']
                m0_3=params['m_min']
                
        k1=1./(piece1+piece2*k2_over_k1+piece3*k3_over_k2*k2_over_k1)#sample size normalized to 1
        k2=k1*k2_over_k1
        k3=k2*k3_over_k2

        #get fraction of sample within each piece
        f1=k1*piece1
        f2=k2*piece2
        f3=k3*piece3

        for i in range(0,len(ran1)):
                
            if ran1[i]<f1:
                mass.append((params['m_min']**(1.-params['alpha1'])+ran1[i]*(1.-params['alpha1'])/k1)**(1./(1.-params['alpha1'])))
            elif ((ran1[i]>=f1)&(ran1[i]<(f1+f2))):
                mass.append((m0_2**(1.-params['alpha2'])+(1.-params['alpha2'])/k2*(ran1[i]-f1))**(1./(1.-params['alpha2'])))
            elif ran1[i]>=(f1+f2):
                mass.append((m0_3**(1.-params['alpha3'])+(1.-params['alpha3'])/k3*(ran1[i]-f1-f2))**(1./(1.-params['alpha3'])))
            else:
                raise ValueError('something wrong in sampling Kroupa IMF')
                
        mass=np.array(mass)

        def kroupa_func(x):
            if ((type(x) is list)|(type(x) is np.ndarray)):
                val=[]
                for i in range(0,len(x)):
                    if x[i]<params['m1_break']:
                        val.append(k1*x[i]**-params['alpha1'])
                    elif ((x[i]>=params['m1_break'])&(x[i]<params['m2_break'])):
                        val.append(k2*x[i]**-params['alpha2'])
                    elif x[i]>=params['m2_break']:
                        val.append(k3*x[i]**-params['alpha3'])
                    else:
                        raise ValueError('problem in kroupa_func')
                val=np.array(val)
                
            elif ((type(x) is float)|(type(x) is int)):
                if x<params['m1_break']:
                    val=k1*x**-params['alpha1']
                elif ((x>=params['m1_break'])&(x<params['m2_break'])):
                    val=k2*x**-params['alpha2']
                elif x>=params['m2_break']:
                    val=k3*x**-params['alpha3']
                else:
                    raise ValueError('problem in kroupa_func')
            else:
                raise TypeError('type error in kroupa func')
                
            return val
        
        return imf(model=model,mass=mass,alpha1=params['alpha1'],alpha2=params['alpha2'],alpha3=params['alpha3'],m1_break=params['m1_break'],m2_break=params['m2_break'],m_min=params['m_min'],m_max=params['m_max'],k1=k1,k2=k2,k3=k3,func=kroupa_func)

    if model=='bpl':#sample from broken power law, 2 separate power laws with indices -alpha1, -alpha2, break mass at m_break

        if not 'alpha1' in params:
            params['alpha1']=1.3
        if not 'alpha2' in params:
            params['alpha2']=2.3
        if not 'm_break' in params:
            params['m_break']=0.5
            
        mass=[]

        #get normalization constant for each of three pieces
        k2_over_k1=params['m_break']**(params['alpha2']-params['alpha1'])
        
        if params['m_min']<params['m_break']:
            
            if params['m_break']<=params['m_max']:
                
                piece1=(params['m_break']**(1.-params['alpha1'])-params['m_min']**(1.-params['alpha1']))/(1.-params['alpha1'])
                piece2=(params['m_max']**(1.-params['alpha2'])-params['m_break']**(1.-params['alpha2']))/(1.-params['alpha2'])

                m0_2=params['m_break']
                                
            if params['m_break']>params['m_max']:
                
                piece1=(params['m_max']**(1.-params['alpha1'])-params['m_min']**(1.-params['alpha1']))/(1.-params['alpha1'])
                piece2=0.
                
                m0_2=params['m_break']
                
        if params['m_break']<=params['m_min']:
            
            if params['m_break']<=params['m_max']:
                piece1=0.
                piece2=(params['m_max']**(1.-params['alpha2'])-params['m_min']**(1.-params['alpha2']))/(1.-params['alpha2'])
                
                m0_2=params['m_min']
                
        k1=1./(piece1+piece2*k2_over_k1)#sample size normalized to 1
        k2=k1*k2_over_k1

        #get fraction of sample within each piece
        f1=k1*piece1
        f2=k2*piece2

        for i in range(0,len(ran1)):
                
            if ran1[i]<f1:
                mass.append((params['m_min']**(1.-params['alpha1'])+ran1[i]*(1.-params['alpha1'])/k1)**(1./(1.-params['alpha1'])))
            elif ((ran1[i]>=f1)&(ran1[i]<(f1+f2))):
                mass.append((m0_2**(1.-params['alpha2'])+(1.-params['alpha2'])/k2*(ran1[i]-f1))**(1./(1.-params['alpha2'])))
            else:
                raise ValueError('something wrong in sampling BPL IMF')
                
        mass=np.array(mass)

        def bpl_func(x):
            if ((type(x) is list)|(type(x) is np.ndarray)):
                val=[]
                for i in range(0,len(x)):
                    if x[i]<params['m_break']:
                        val.append(k1*x[i]**-params['alpha1'])
                    elif x[i]>=params['m_break']:
                        val.append(k2*x[i]**-params['alpha2'])
                    else:
                        raise ValueError('problem in bpl_func')
                val=np.array(val)
                
            elif ((type(x) is float)|(type(x) is int)):
                if x<params['m_break']:
                    val=k1*x**-params['alpha1']
                elif x>=params['m_break']:
                    val=k2*x**-params['alpha2']
                else:
                    raise ValueError('problem in bpl_func')
            else:
                raise TypeError('type error in bpl func')
                
            return val
        
        return imf(model=model,mass=mass,alpha1=params['alpha1'],alpha2=params['alpha2'],m_break=params['m_break'],m_min=params['m_min'],m_max=params['m_max'],k1=k1,k2=k2,func=bpl_func)

def sample_orbit_2body(f_period,**params):#f_period is time of observation / period, with f_period=0 at pericenter.  Can handle f_period < 1e10.

    class orbit_2body:
        
        def __init__(self,semimajor_axis=None,eccentricity=None,mass_primary=None,mass_secondary=None,energy=None,angular_momentum=None,f_period=None,time=None,period=None,eta=None,theta=None,r_xyz=None,r_sph=None,v_xyz=None,v_sph=None,r1_xyz=None,v1_xyz=None,r1_sph=None,v1_sph=None,r2_xyz=None,v2_xyz=None,r2_sph=None,v2_sph=None,inclination=None,longitude=None,r_obs_xyz=None,v_obs_xyz=None,r_obs_sph=None,v_obs_sph=None,r1_obs_xyz=None,v1_obs_xyz=None,r1_obs_sph=None,v1_obs_sph=None,r2_obs_xyz=None,v2_obs_xyz=None,r2_obs_sph=None,v2_obs_sph=None,rot_matrix=None):
            self.semimajor_axis=semimajor_axis #AU
            self.eccentricity=eccentricity
            self.mass_primary=mass_primary #Msun
            self.mass_secondary=mass_secondary #Msun
            self.energy=energy #total orbital energy per (reduced) mass, units of AU^2 / yr^2
            self.angular_momentum=angular_momentum #total orbital angular momentum per (reduced) mass, units of AU^2/yr
            self.f_period=f_period #time/period
            self.time=time #time sinze t=0 at theta=0, yr
            self.period=period #orbital period, yr
            self.eta=eta #eccentric anomaly (radians)
            self.theta=theta #true anomaly (radians)
            self.r_xyz=r_xyz #reduced mass position in CM frame, AU
            self.v_xyz=v_xyz #reduced mass velocity in CM frame, AU/yr
            self.r_sph=r_sph #reduced mass position (r,longitude,inclination), in (AU, radians, radians); longitude is azimuthal angle 
            self.v_sph=v_sph #reduced mass velocity (v_r,v_longidue,v_inclination) in AU/yr
            self.r1_xyz=r1_xyz #particle 1 position, AU
            self.v1_xyz=v1_xyz #particle 1 velocity, AU/yr
            self.r2_xyz=r2_xyz #particle 2 position, AU
            self.v2_xyz=v2_xyz #particle 2 velocity, AU/yr
            self.inclination=inclination #inclination defined by observer's position, radians
            self.longitude=longitude #azimuthal angle defined by observer's position, radians
            self.r_obs_xyz=r_obs_xyz
            self.v_obs_xyz=v_obs_xyz
            self.r1_obs_xyz=r1_obs_xyz            
            self.v1_obs_xyz=v1_obs_xyz
            self.r2_obs_xyz=r2_obs_xyz
            self.v2_obs_xyz=v2_obs_xyz
            self.rot_matrix=rot_matrix

    #default is Sun/Earth orbit
    if not 'period' in params:
        params['period']=1.*u.yr
    if not 'eccentricity' in params:
        params['eccentricity']=0.
    if not 'mass_primary' in params:
        params['mass_primary']=1.*u.M_sun
    if not 'mass_ratio' in params:
        params['mass_ratio']=1.

    #if any of f_period, mass_primary, mass_secondary, period, eccentricity, inclination, longitude are input as scalars, convert to arrays of same length as f_period (if f_period is input as scalar, first make it array of length 1)
    
    f_period=np.array(f_period)
    params['eccentricity']=np.array(params['eccentricity'])
    params['inclination']=np.array(params['inclination'])*params['inclination'].unit
    params['longitude']=np.array(params['longitude'])*params['longitude'].unit
    params['mass_primary']=np.array(params['mass_primary'])*params['mass_primary'].unit
    params['mass_ratio']=np.array(params['mass_ratio'])
    params['period']=np.array(params['period'])*params['period'].unit

    if np.size(f_period)==1:
        f_period=np.array([f_period]).reshape(1)
    if np.size(params['eccentricity'])==1:
        params['eccentricity']=np.full(len(f_period),np.array(params['eccentricity']).reshape(1))
    if np.size(params['inclination'])==1:
        params['inclination']=np.full(len(f_period),np.array(params['inclination']).reshape(1))*params['inclination'].unit
    if np.size(params['longitude'])==1:
        params['longitude']=np.full(len(f_period),np.array(params['longitude']).reshape(1))*params['longitude'].unit
    if np.size(params['mass_primary'])==1:
        params['mass_primary']=np.full(len(f_period),np.array(params['mass_primary']).reshape(1))*params['mass_primary'].unit
    if np.size(params['mass_ratio'])==1:
        params['mass_ratio']=np.full(len(f_period),np.array(params['mass_ratio']).reshape(1))
    if np.size(params['period'])==1:
        params['period']=np.full(len(f_period),np.array(params['period']).reshape(1))*params['period'].unit

    if not ((len(params['eccentricity'])==len(f_period))&(len(params['inclination'])==len(f_period))&(len(params['longitude'])==len(f_period))&(len(params['mass_primary'])==len(f_period))&(len(params['mass_ratio'])==len(f_period))&(len(params['period'])==len(f_period))):
        raise ValueError("if input as lists or arrays with size>1, 'eccentricity', 'inclination', 'longitude', 'mass_primary', 'mass_ratio', 'period' must all be of same length.  If any are input as lists with one element, will be understood to apply to all times")
            
    g=4*np.pi**2*u.AU**3/u.yr**2/u.M_sun#keplerian units (AU, yr, Msun)

    mass_secondary=params['mass_primary']*params['mass_ratio']
    mass=params['mass_primary']+mass_secondary
    semimajor_axis=(params['period']**2*g*mass/4.*np.pi**2)**(1./3.) #AU, assuming period in yr and mass in Msun
    #period=np.sqrt(params['semimajor_axis']**3/(params['mass_primary']+mass_secondary)) #yr, assuming semimajor axis in AU and mass in Msun
    energy=-g*mass/2/semimajor_axis
    angular_momentum=np.sqrt(g*mass*(1.-params['eccentricity']**2))
    #function to solve for eta (eccentric anomaly, defined in Ch. 3.1 of Binney/Tremaine 2008) as a function of f_period=t/period
    def find_eta(x,eccentricity,f_period):
        return (x-eccentricity*np.sin(x))/2/np.pi-f_period

    #compute eta (eccentric anomaly) according to f_period = t/period array
    low=0.
    high=2.*np.pi
    
    f_period_eff=np.zeros(len(f_period ))
    for i in range(0,len(f_period)):
        f_period_eff[i]=np.modf(f_period[i])[0] #fraction of period after removing all completed periods (so eta remains within interval (0, 2pi))
        
    eta=np.zeros(len(f_period))
    for i in range(0,len(f_period)):
        eta[i]=scipy.optimize.brentq(find_eta,low,high,args=(params['eccentricity'][i],f_period_eff[i]))#this is eccentric anomaly
    
    #use Eq. 3.28a from Binney/Tremaine 2008 to calculate r as function of eta.  vector(r) = vector(r2) - vector(r1) represents separation between particles 1 and 2.
    r=semimajor_axis*(1.-params['eccentricity']*np.cos(eta))# has same units as semi-major axis
    #use Eq. 3.326 from Binney/Tremaine 2008 to convert eccentric anomaly into true anomaly theta.
    theta=np.arccos((np.cos(eta)-params['eccentricity'])/(1.-params['eccentricity']*np.cos(eta)))
    #if ((type(eta)==float)|(type(eta)==np.float64)|(type(eta)==np.float32)):
        #if eta>=np.pi:
            #theta=2.*np.pi-theta#fudge for quadrant problem
    #else:
    theta[eta>=np.pi]=2.*np.pi-theta[eta>=np.pi]#fudge for quadrant problem

    #get velocity of reduced mass, AU/yr
          
    v=2*np.pi*semimajor_axis/params['period']*np.sqrt(2.*(1.+params['eccentricity']*np.cos(theta))/(1.-params['eccentricity']**2)-1.)
    vr=2.*np.pi*semimajor_axis/params['period']*params['eccentricity']*np.sin(theta)/np.sqrt(1.-params['eccentricity']**2)
    vtheta=2.*np.pi*semimajor_axis/params['period']*(1.+params['eccentricity']*np.cos(theta))/np.sqrt(1.-params['eccentricity']**2)
                   
    #get Cartesian coordinates of separation vector and velocity of reduced mass
    x=r*np.cos(theta)
    y=r*np.sin(theta)
    vx=vr*np.cos(theta)-vtheta*np.sin(theta)
    vy=vr*np.sin(theta)+vtheta*np.cos(theta)
    
    #transform r and v into position and velocity vectors (and x,y components) for real particles 1 and 2
    if len(np.where(params['mass_ratio']>1)[0])>0:
        raise ValueError('must have mass_primary > mass_secondary')
    
    trans1,trans2=-params['mass_ratio']/(1.+params['mass_ratio']),1./(1.+params['mass_ratio'])
    
    r1,r2=trans1*r,trans2*r #AU
    v1,v2=trans1*v,trans2*v #AU/yr
                   
    x1,y1=trans1*x,trans1*y #AU
    vx1,vy1=trans1*vx,trans1*vy #AU/yr
                   
    x2,y2=trans2*x,trans2*y #AU
    vx2,vy2=trans2*vx,trans2*vy #AU/yr

    z,vz=x-x,x-x #orbit is confined to xy plane

    r_xyz=np.array((x,y,z)).T*r.unit
    r_sph=np.array((r,theta,z)).T
    v_xyz=np.array((vx,vy,z)).T*v.unit
    v_sph=np.array((vr,vtheta,vz)).T
    r1_xyz=np.array((x1,y1,z)).T*r.unit
    v1_xyz=np.array((vx1,vy1,z)).T*v.unit
    r2_xyz=np.array((x2,y2,z)).T*r.unit
    v2_xyz=np.array((vx2,vy2,vz)).T*v.unit

    r_obs_xyz=np.zeros(np.shape(r_xyz))*r.unit
    r1_obs_xyz=np.zeros(np.shape(r_xyz))*r.unit
    r2_obs_xyz=np.zeros(np.shape(r_xyz))*r.unit
    v_obs_xyz=np.zeros(np.shape(r_xyz))*v.unit
    v1_obs_xyz=np.zeros(np.shape(r_xyz))*v.unit
    v2_obs_xyz=np.zeros(np.shape(r_xyz))*v.unit
    
    if 'inclination' in params:

        rot_matrix=[]
        
        rot_alpha=params['longitude'].to(u.rad).value #rotation about z axis, in direction of arc from +x to +y (radians)
        rot_beta=0. #rotation about y axis, in direction of arc from +z to +x (radians)
        rot_gamma=params['inclination'].to(u.rad).value #rotation about x axis, in direction of arc from +y to +z (radians)

        if ((len(params['longitude'])!=len(f_period))|(len(params['inclination'])!=len(f_period))):
            raise ValueError("'inclination' and 'longitude' must have same length as 'f_period' array")
        for i in range(0,len(params['longitude'])):
            rot_matrix.append(get_rot_matrix(rot_alpha[i],rot_beta,rot_gamma[i]))

        for i in range(0,len(rot_matrix)):
            r_obs_xyz.value[i]=rot_matrix[i].apply(r_xyz.value[i])
            r1_obs_xyz.value[i]=rot_matrix[i].apply(r1_xyz.value[i])
            r2_obs_xyz.value[i]=rot_matrix[i].apply(r2_xyz.value[i])    
            v_obs_xyz.value[i]=rot_matrix[i].apply(v_xyz.value[i])
            v1_obs_xyz.value[i]=rot_matrix[i].apply(v1_xyz.value[i])
            v2_obs_xyz.value[i]=rot_matrix[i].apply(v2_xyz.value[i])

    return orbit_2body(semimajor_axis=semimajor_axis,eccentricity=params['eccentricity'],mass_primary=params['mass_primary'],mass_secondary=mass_secondary,energy=energy,angular_momentum=angular_momentum,inclination=params['inclination'],longitude=params['longitude'],f_period=f_period,time=f_period*params['period'],period=params['period'],eta=eta,theta=theta,r_xyz=r_xyz,r_sph=r_sph,v_xyz=v_xyz,v_sph=v_sph,r1_xyz=r1_xyz,v1_xyz=v1_xyz,r2_xyz=r2_xyz,v2_xyz=v2_xyz,r_obs_xyz=r_obs_xyz,r1_obs_xyz=r1_obs_xyz,r2_obs_xyz=r2_obs_xyz,v_obs_xyz=v_obs_xyz,v1_obs_xyz=v1_obs_xyz,v2_obs_xyz=v2_obs_xyz,rot_matrix=rot_matrix)

def sample_normal_truncated(**params):
    if not 'size' in params:
        params['size']=1
    if not 'min_value' in params:
        params['min_value']=-np.inf
    if not 'max_value' in params:
        params['max_value']=np.inf
    if not 'loc' in params:
        params['loc']=0.
    if not 'scale' in params:
        params['scale']=1.

    return params['loc']-np.sqrt(2.*params['scale']**2)*scipy.special.erfinv(scipy.special.erf((params['loc']-params['min_value'])/np.sqrt(2.*params['scale']**2))-np.random.uniform(size=params['size'],low=0.,high=1.)*(scipy.special.erf((params['loc']-params['min_value'])/np.sqrt(2.*params['scale']**2))-scipy.special.erf((params['loc']-params['max_value'])/np.sqrt(2.*params['scale']**2))))

def sample_thermal(**params):
    if not 'size' in params:
        params['size']=1
    return np.sqrt(np.random.uniform(size=params['size'],low=0.,high=1.))

def sample_inclination(**params):
    if not 'size' in params:
        params['size']=1
    ran1=np.random.uniform(size=params['size'],low=0.,high=1.)
    ran2=np.random.uniform(size=params['size'],low=0.,high=1.)
    inclination=np.arccos((1.-2*ran1))
    change=np.where(ran2>0.5)[0]
    inclination[change]=inclination[change]+np.pi
    return inclination

def get_rot_matrix(alpha,beta,gamma):
    #alpha is rotation about z axis, from +x to +y ('yaw' in radians)
    #beta is rotation about y axis, from +z to +x ('pitch' in radians)
    #gamma is rotation about x axis, from +y to +z ('roll' in radians)
    #rotations are performed in the order gamma, beta, alpha
    r11=np.cos(alpha)*np.cos(beta)
    r12=np.cos(alpha)*np.sin(beta)*np.sin(gamma)-np.sin(alpha)*np.cos(gamma)
    r13=np.cos(alpha)*np.sin(beta)*np.cos(gamma)+np.sin(alpha)*np.sin(gamma)
    r21=np.sin(alpha)*np.cos(beta)
    r22=np.sin(alpha)*np.sin(beta)*np.sin(gamma)+np.cos(alpha)*np.cos(gamma)
    r23=np.sin(alpha)*np.sin(beta)*np.cos(gamma)-np.cos(alpha)*np.sin(gamma)
    r31=-np.sin(beta)
    r32=np.cos(beta)*np.sin(gamma)
    r33=np.cos(beta)*np.cos(gamma)
    return Rotation.from_matrix([[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]])

def sample_combine(sample_r2d,sample_imf,sample_binary,sample_orbit):

    class sample_final:
        
        def __init__(self,r_xyz=None,mass=None,item=None,companion=None):
            self.r_xyz=r_xyz #AU
            self.mass=mass
            self.item=item
            self.companion=companion

    r_xyz,mass,item,companion=[],[],[],[]
    j=0
    for i in range(0,len(sample_r2d.r_xyz)):
        if sample_binary[i]:
            if sample_imf.mass[i] != sample_orbit.mass_primary[j]:
                raise ValueError ("problem with sample_final masses")
            
            r_xyz.append((sample_r2d.r_xyz[i]+sample_orbit.r1_obs_xyz[j]).tolist())
            mass.append(sample_orbit.mass_primary[j])
            item.append('primary')
            companion.append(i+j+1)
            
            r_xyz.append((sample_r2d.r_xyz[i]+sample_orbit.r2_obs_xyz[j]).tolist())
            mass.append(sample_orbit.mass_secondary[j])
            item.append('secondary')
            companion.append(i+j)

            j+=1
            
        else:
            
            r_xyz.append(sample_r2d.r_xyz[i].tolist())
            mass.append(sample_imf.mass[i])
            item.append('single')
            companion.append(-999)
            
    return sample_final(r_xyz=np.array(r_xyz),mass=np.array(mass),item=np.array(item),companion=np.array(companion,dtype=int))

def animation_2body_r(sample_orbit,animation_filename):
    class animate_r(object):#based on https://stackoverflow.com/questions/9401658/how-to-animate-a-scatter-plot

        def __init__(self):
        
            self.stream=self.data_stream()

            #Set up figure and axes.
            self.fig=plt.figure(1)

            self.ax1=self.fig.add_subplot(221,aspect='equal')
            self.ax2=self.fig.add_subplot(222,aspect='equal')
            self.ax3=self.fig.add_subplot(223,aspect='equal')
            self.ax4=self.fig.add_subplot(224,aspect='equal')
            self.fig.subplots_adjust(left=-0,wspace=-0.3,hspace=0.25)

            #set up funcanimation
            if animation_filename==None:
                save_animation_to_file=False
            else:
                save_animation_to_file=True
            if save_animation_to_file:
                self.ani = animation.FuncAnimation(self.fig, self.update, interval=5,init_func=self.setup_plot, blit=False,save_count=1*len(sample_orbit.r_xyz))
                #writer=animation.PillowWriter(fps=60)
                #f='2body_r.gif'
                writer=animation.FFMpegWriter(fps=60)
                f=animation_filename+'.mp4'
                self.ani.save(f,writer=writer)
            else:
                self.ani = animation.FuncAnimation(self.fig, self.update, interval=5,init_func=self.setup_plot, blit=False,save_count=1*len(sample_orbit.r_xyz))

        def setup_plot(self):
            
            first=True
            x,y,z,x1,y1,z1,x2,y2,z2,x_obs,y_obs,z_obs,x1_obs,y1_obs,z1_obs,x2_obs,y2_obs,z2_obs,t= next(self.stream).T

            lim=np.max(np.abs(np.concatenate([sample_orbit.r_xyz.T[0].value,sample_orbit.r_xyz.T[1].value,sample_orbit.r1_xyz.T[0].value,sample_orbit.r1_xyz.T[1].value,sample_orbit.r2_xyz.T[0].value,sample_orbit.r2_xyz.T[1].value])))*1.2/sample_orbit.semimajor_axis[0].value

            self.ax1.plot(sample_orbit.r_xyz.T[0].value/sample_orbit.semimajor_axis.value,sample_orbit.r_xyz.T[1].value/sample_orbit.semimajor_axis.value,color='k',lw=1)
            self.ax1.plot(sample_orbit.r1_xyz.T[0].value/sample_orbit.semimajor_axis.value,sample_orbit.r1_xyz.T[1].value/sample_orbit.semimajor_axis.value,color='b',lw=1)
            self.ax1.plot(sample_orbit.r2_xyz.T[0].value/sample_orbit.semimajor_axis.value,sample_orbit.r2_xyz.T[1].value/sample_orbit.semimajor_axis.value,color='r',lw=1)
            self.ax2.plot(sample_orbit.r_xyz.T[1].value/sample_orbit.semimajor_axis.value,sample_orbit.r_xyz.T[2].value/sample_orbit.semimajor_axis.value,color='k',lw=1)
            self.ax2.plot(sample_orbit.r1_xyz.T[1].value/sample_orbit.semimajor_axis.value,sample_orbit.r1_xyz.T[2].value/sample_orbit.semimajor_axis.value,color='b',lw=1)
            self.ax2.plot(sample_orbit.r2_xyz.T[1].value/sample_orbit.semimajor_axis.value,sample_orbit.r2_xyz.T[2].value/sample_orbit.semimajor_axis.value,color='r',lw=1)
            self.ax3.plot(sample_orbit.r_obs_xyz.T[0].value/sample_orbit.semimajor_axis.value,sample_orbit.r_obs_xyz.T[1].value/sample_orbit.semimajor_axis.value,color='k',lw=1)
            self.ax3.plot(sample_orbit.r1_obs_xyz.T[0].value/sample_orbit.semimajor_axis.value,sample_orbit.r1_obs_xyz.T[1].value/sample_orbit.semimajor_axis.value,color='b',lw=1)
            self.ax3.plot(sample_orbit.r2_obs_xyz.T[0].value/sample_orbit.semimajor_axis.value,sample_orbit.r2_obs_xyz.T[1].value/sample_orbit.semimajor_axis.value,color='r',lw=1)
            self.ax4.plot(sample_orbit.r_obs_xyz.T[1].value/sample_orbit.semimajor_axis.value,sample_orbit.r_obs_xyz.T[2].value/sample_orbit.semimajor_axis.value,color='k',lw=1)
            self.ax4.plot(sample_orbit.r1_obs_xyz.T[1].value/sample_orbit.semimajor_axis.value,sample_orbit.r1_obs_xyz.T[2].value/sample_orbit.semimajor_axis.value,color='b',lw=1)
            self.ax4.plot(sample_orbit.r2_obs_xyz.T[1].value/sample_orbit.semimajor_axis.value,sample_orbit.r2_obs_xyz.T[2].value/sample_orbit.semimajor_axis.value,color='r',lw=1)
        
            self.ax1.text(-lim*0.9,lim*0.88,r'$m_2/m_1=$'+str("{:.2f}".format(round(sample_orbit.mass_secondary[0].value/sample_orbit.mass_primary[0].value,3))),fontsize=5) #after rotation
            self.ax1.text(-lim*0.9,lim*0.76,r'$e=$'+str("{:.2f}".format(round(sample_orbit.eccentricity[0],3))),fontsize=5)
            self.ax3.text(-lim*0.9,lim*0.88,r'$i=$'+str("{:.2f}".format(round(sample_orbit.inclination[0].to(u.rad).value*180/np.pi,2)))+r'$^{\circ}$',fontsize=5)
            self.ax3.text(-lim*0.9,lim*0.76,r'$l=$'+str("{:.2f}".format(round(sample_orbit.longitude[0].to(u.rad).value*180/np.pi,2)))+r'$^{\circ}$',fontsize=5)

            self.ax1.set_xlabel(r'$x/a$',fontsize=7)
            self.ax1.set_ylabel(r'$y/a$',fontsize=7,labelpad=-3)
            
            self.ax2.set_xlabel(r'$y/a$',fontsize=7)
            self.ax2.set_ylabel(r'$z/a$',fontsize=7,labelpad=-3)

            self.ax3.set_xlabel(r'$x_{\rm obs}/a$',fontsize=7)
            self.ax3.set_ylabel(r'$y_{\rm obs}/a$',fontsize=7,labelpad=-3)

            self.ax4.set_xlabel(r'$y_{\rm obs}/a$',fontsize=7)
            self.ax4.set_ylabel(r'$z_{\rm obs}/a$',fontsize=7,labelpad=-3)
        
            for ax in [self.ax1,self.ax2,self.ax3,self.ax4]:
                ax.axvline(0,ls=':',color='k')
                ax.axhline(0,ls=':',color='k')
                ax.set_xlim([-lim,lim])
                ax.set_ylim([-lim,lim])
                ax.tick_params(labelsize=7)

            if first:
                self.ax1.scatter([-np.inf],[-np.inf],s=15,edgecolor='k',facecolor='none',label='reduced mass')
                self.ax1.scatter([-np.inf],[-np.inf],s=20,color='b',label='particle 1')
                self.ax1.scatter([-np.inf],[-np.inf],s=10,color='r',label='particle 2')
                self.ax1.legend(fontsize=5)
                first=False
                
            #return the artists 
            self.xy_scat = self.ax1.scatter(x, y, edgecolor='k', facecolor='none',s=15)
            self.xy1_scat = self.ax1.scatter(x1, y1, c='b', s=20)
            self.xy2_scat = self.ax1.scatter(x2, y2, c='r', s=10)
            self.yz_scat = self.ax2.scatter(y, z, edgecolor='k', facecolor='none',s=15)
            self.yz1_scat = self.ax2.scatter(y1, z1, c='b', s=20)
            self.yz2_scat = self.ax2.scatter(y2, z2, c='r', s=10)
            self.xy_obs_scat = self.ax3.scatter(x_obs, y_obs, edgecolor='k', facecolor='none',s=15)
            self.xy1_obs_scat = self.ax3.scatter(x1_obs, y1_obs, c='b', s=20)
            self.xy2_obs_scat = self.ax3.scatter(x2_obs, y2_obs, c='r', s=10)
            self.yz_obs_scat = self.ax4.scatter(y_obs, z_obs, edgecolor='k', facecolor='none',s=15)
            self.yz1_obs_scat = self.ax4.scatter(y1_obs, z1_obs, c='b', s=20)
            self.yz2_obs_scat = self.ax4.scatter(y2_obs, z2_obs, c='r', s=10)
            self.t = self.ax1.text(-lim*0.9,lim*0.64,r'time / period $=$'+str("{:.2f}".format(round(t[0],2))),fontsize=5)
            
            return self.xy_scat,self.xy1_scat,self.xy2_scat,self.xy_obs_scat,self.xy1_obs_scat,self.xy2_obs_scat,self.yz_scat,self.yz1_scat,self.yz2_scat,self.yz_obs_scat,self.yz1_obs_scat,self.yz2_obs_scat,self.t
    
        def data_stream(self):
            
            i=0
            while True:
                x=sample_orbit.r_xyz.T[0][i].value/sample_orbit.semimajor_axis[i].value
                y=sample_orbit.r_xyz.T[1][i].value/sample_orbit.semimajor_axis[i].value
                z=sample_orbit.r_xyz.T[2][i].value/sample_orbit.semimajor_axis[i].value
                x1=sample_orbit.r1_xyz.T[0][i].value/sample_orbit.semimajor_axis[i].value
                y1=sample_orbit.r1_xyz.T[1][i].value/sample_orbit.semimajor_axis[i].value
                z1=sample_orbit.r1_xyz.T[2][i].value/sample_orbit.semimajor_axis[i].value
                x2=sample_orbit.r2_xyz.T[0][i].value/sample_orbit.semimajor_axis[i].value
                y2=sample_orbit.r2_xyz.T[1][i].value/sample_orbit.semimajor_axis[i].value
                z2=sample_orbit.r2_xyz.T[2][i].value/sample_orbit.semimajor_axis[i].value
                x_obs=sample_orbit.r_obs_xyz.T[0][i].value/sample_orbit.semimajor_axis[i].value
                y_obs=sample_orbit.r_obs_xyz.T[1][i].value/sample_orbit.semimajor_axis[i].value
                z_obs=sample_orbit.r_obs_xyz.T[2][i].value/sample_orbit.semimajor_axis[i].value
                x1_obs=sample_orbit.r1_obs_xyz.T[0][i].value/sample_orbit.semimajor_axis[i].value
                y1_obs=sample_orbit.r1_obs_xyz.T[1][i].value/sample_orbit.semimajor_axis[i].value
                z1_obs=sample_orbit.r1_obs_xyz.T[2][i].value/sample_orbit.semimajor_axis[i].value
                x2_obs=sample_orbit.r2_obs_xyz.T[0][i].value/sample_orbit.semimajor_axis[i].value
                y2_obs=sample_orbit.r2_obs_xyz.T[1][i].value/sample_orbit.semimajor_axis[i].value
                z2_obs=sample_orbit.r2_obs_xyz.T[2][i].value/sample_orbit.semimajor_axis[i].value
                t=sample_orbit.time[i].value/sample_orbit.period[i].value
                i+=1
                if i==len(sample_orbit.r_xyz):#infinite loop 
                    i=0
                yield np.c_[x,y,z,x1,y1,z1,x2,y2,z2,x_obs,y_obs,z_obs,x1_obs,y1_obs,z1_obs,x2_obs,y2_obs,z2_obs,t]

        def update(self, i):
        
            data = next(self.stream)

            self.xy_scat.set_offsets([[data[0][0],data[0][1]]])
            self.yz_scat.set_offsets([[data[0][1],data[0][2]]])
            self.xy1_scat.set_offsets([[data[0][3],data[0][4]]])
            self.yz1_scat.set_offsets([[data[0][4],data[0][5]]])
            self.xy2_scat.set_offsets([[data[0][6],data[0][7]]])
            self.yz2_scat.set_offsets([[data[0][7],data[0][8]]])

            self.xy_obs_scat.set_offsets([[data[0][9],data[0][10]]])
            self.yz_obs_scat.set_offsets([[data[0][10],data[0][11]]])
            self.xy1_obs_scat.set_offsets([[data[0][12],data[0][13]]])
            self.yz1_obs_scat.set_offsets([[data[0][13],data[0][14]]])
            self.xy2_obs_scat.set_offsets([[data[0][15],data[0][16]]])
            self.yz2_obs_scat.set_offsets([[data[0][16],data[0][17]]])
            
            self.t.set_text('time / period='+str("{:.2f}".format(round(data[:,18:19][0][0],2))))

            #return updated artist
            return self.xy_scat,self.xy1_scat,self.xy2_scat,self.xy_obs_scat,self.xy1_obs_scat,self.xy2_obs_scat,self.yz_scat,self.yz1_scat,self.yz2_scat,self.yz_obs_scat,self.yz1_obs_scat,self.yz2_obs_scat,self.t
    a=animate_r()
    plt.show()



def animation_2body_v(sample_orbit,animation_filename):

    class animate_v(object):#based on https://stackoverflow.com/questions/9401658/how-to-animate-a-scatter-plot

        def __init__(self):
            self.stream = self.data_stream()

            self.fig=plt.figure(1)
            self.ax1=self.fig.add_subplot(231)
            self.ax2=self.fig.add_subplot(232)
            self.ax3=self.fig.add_subplot(233)
            self.ax4=self.fig.add_subplot(234)
            self.ax5=self.fig.add_subplot(235)
            self.ax6=self.fig.add_subplot(236)
            self.fig.subplots_adjust(wspace=0.45,hspace=0.25)

            if animation_filename==None:
                save_animation_to_file=False
            else:
                save_animation_to_file=True
                
            if save_animation_to_file:
                self.ani = animation.FuncAnimation(self.fig, self.update, interval=5,init_func=self.setup_plot, blit=False,save_count=10*len(sample_orbit.r_xyz))
                writer=animation.FFMpegWriter(fps=60)
                f=aniation_filename+'.mp4'
                self.ani.save(f,writer=writer)
            else:
                self.ani = animation.FuncAnimation(self.fig, self.update, interval=5,init_func=self.setup_plot, blit=False,save_count=1*len(sample_orbit.r_xyz))

        def setup_plot(self):

            first=True
            t,v_x,v1_x,v2_x,v_y,v1_y,v2_y,v_z,v1_z,v2_z,v_obs_x,v1_obs_x,v2_obs_x,v_obs_y,v1_obs_y,v2_obs_y,v_obs_z,v1_obs_z,v2_obs_z,t_text= next(self.stream).T

            lim=np.max(np.abs(np.concatenate([sample_orbit.v_xyz.T[0].value,sample_orbit.v_xyz.T[1].value,sample_orbit.v_xyz.T[2].value,sample_orbit.v1_xyz.T[0].value,sample_orbit.v1_xyz.T[1].value,sample_orbit.v1_xyz.T[2].value,sample_orbit.v2_xyz.T[0].value,sample_orbit.v2_xyz.T[1].value,sample_orbit.v2_xyz.T[2].value])))*1.1/np.sqrt(-2.*sample_orbit.energy[0].value)

            if first:
                self.ax3.scatter([-np.inf],[-np.inf],s=15,edgecolor='k',facecolor='none',label='reduced mass')
                self.ax3.scatter([-np.inf],[-np.inf],s=20,color='b',label='particle 1')
                self.ax3.scatter([-np.inf],[-np.inf],s=10,color='r',label='particle 2')
                first=False
        
            self.ax1.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v_xyz.T[0].value/np.sqrt(-2.*sample_orbit.energy.value),color='k',lw=1)
            self.ax1.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v1_xyz.T[0].value/np.sqrt(-2.*sample_orbit.energy.value),color='b',lw=1)
            self.ax1.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v2_xyz.T[0].value/np.sqrt(-2.*sample_orbit.energy.value),color='r',lw=1)
        
            self.ax2.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v_xyz.T[1].value/np.sqrt(-2.*sample_orbit.energy.value),color='k',lw=1)
            self.ax2.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v1_xyz.T[1].value/np.sqrt(-2.*sample_orbit.energy.value),color='b',lw=1)
            self.ax2.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v2_xyz.T[1].value/np.sqrt(-2.*sample_orbit.energy.value),color='r',lw=1)
        
            self.ax3.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v_xyz.T[2].value/np.sqrt(-2.*sample_orbit.energy.value),color='k',lw=1)
            self.ax3.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v1_xyz.T[2].value/np.sqrt(-2.*sample_orbit.energy.value),color='b',lw=1)
            self.ax3.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v2_xyz.T[2].value/np.sqrt(-2.*sample_orbit.energy.value),color='r',lw=1)

            self.ax4.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v_obs_xyz.T[0].value/np.sqrt(-2.*sample_orbit.energy.value),color='k',lw=1)
            self.ax4.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v1_obs_xyz.T[0].value/np.sqrt(-2.*sample_orbit.energy.value),color='b',lw=1)
            self.ax4.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v2_obs_xyz.T[0].value/np.sqrt(-2.*sample_orbit.energy.value),color='r',lw=1)

            self.ax5.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v_obs_xyz.T[1].value/np.sqrt(-2.*sample_orbit.energy.value),color='k',lw=1)
            self.ax5.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v1_obs_xyz.T[1].value/np.sqrt(-2.*sample_orbit.energy.value),color='b',lw=1)
            self.ax5.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v2_obs_xyz.T[1].value/np.sqrt(-2.*sample_orbit.energy.value),color='r',lw=1)

            self.ax6.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v_obs_xyz.T[2].value/np.sqrt(-2.*sample_orbit.energy.value),color='k',lw=1)
            self.ax6.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v1_obs_xyz.T[2].value/np.sqrt(-2.*sample_orbit.energy.value),color='b',lw=1)
            self.ax6.plot(sample_orbit.time.value/sample_orbit.period.value,sample_orbit.v2_obs_xyz.T[2].value/np.sqrt(-2.*sample_orbit.energy.value),color='r',lw=1)
        
            self.ax1.text(0.05,lim*0.88,r'$m_2/m_1=$'+str("{:.2f}".format(round(sample_orbit.mass_secondary[0].value/sample_orbit.mass_primary[0].value,3))),fontsize=6) #after rotation
            self.ax1.text(0.05,lim*0.76,r'$e=$'+str("{:.2f}".format(round(sample_orbit.eccentricity[0],3))),fontsize=6)
            self.ax4.text(0.05,lim*0.88,r'$i=$'+str("{:.2f}".format(round(sample_orbit.inclination[0].to(u.rad).value*180/np.pi,2)))+r'$^{\circ}$',fontsize=6,horizontalalignment='left')
            self.ax4.text(0.05,lim*0.76,r'$l=$'+str("{:.2f}".format(round(sample_orbit.longitude[0].to(u.rad).value*180/np.pi,2)))+r'$^{\circ}$',fontsize=6,horizontalalignment='left')
            self.ax1.set_ylabel(r'$v_x/\sqrt{-2E}$',fontsize=7)
            self.ax2.set_ylabel(r'$v_y/\sqrt{-2E}$',fontsize=7)
            self.ax3.set_ylabel(r'$v_z/\sqrt{-2E}$',fontsize=7)
            self.ax4.set_ylabel(r'$v_{\rm obs,x}/\sqrt{-2E}$',fontsize=7)
            self.ax5.set_ylabel(r'$v_{\rm obs,y}/\sqrt{-2E}$',fontsize=7)
            self.ax6.set_ylabel(r'$v_{\rm obs,z}/\sqrt{-2E}$',fontsize=7)
        
            for ax in [self.ax1,self.ax2,self.ax3,self.ax4,self.ax5,self.ax6]:
                ax.axvline(0,ls=':',color='k')
                ax.axhline(0,ls=':',color='k')
                ax.set_xlim([0,1])
                ax.set_ylim([-lim,lim])
                ax.set_xlabel('time / period',fontsize=7)
                ax.tick_params(labelsize=7)

            self.v_x_scat=self.ax1.scatter(t, v_x, edgecolor='k', facecolor='none',s=15)
            self.v1_x_scat=self.ax1.scatter(t, v1_x, c='b', s=20)
            self.v2_x_scat=self.ax1.scatter(t, v2_x, c='r', s=10)
            self.v_y_scat=self.ax2.scatter(t, v_y, edgecolor='k', facecolor='none',s=15)
            self.v1_y_scat=self.ax2.scatter(t, v1_y, c='b', s=20)
            self.v2_y_scat=self.ax2.scatter(t, v2_y, c='r', s=10)
            self.v_z_scat=self.ax3.scatter(t, v_z, edgecolor='k',facecolor='none', s=15)
            self.v1_z_scat=self.ax3.scatter(t, v1_z, c='b', s=20)
            self.v2_z_scat=self.ax3.scatter(t, v2_z, c='r', s=10)

            self.v_obs_x_scat=self.ax4.scatter(t, v_obs_x, ec='k',fc='none', s=15)
            self.v1_obs_x_scat=self.ax4.scatter(t, v1_obs_x, c='b', s=20)
            self.v2_obs_x_scat=self.ax4.scatter(t, v2_obs_x, c='r', s=10)
            self.v_obs_y_scat=self.ax5.scatter(t, v_obs_y, ec='k',fc='none', s=15)
            self.v1_obs_y_scat=self.ax5.scatter(t, v1_obs_y, c='b', s=20)
            self.v2_obs_y_scat=self.ax5.scatter(t, v2_obs_y, c='r', s=10)
            self.v_obs_z_scat=self.ax6.scatter(t, v_obs_z, ec='k',fc='none', s=15)
            self.v1_obs_z_scat=self.ax6.scatter(t, v1_obs_z, c='b', s=20)
            self.v2_obs_z_scat=self.ax6.scatter(t, v2_obs_z, c='r', s=10)
        
            self.t_text=self.ax1.text(0.05,lim*0.64,r'time / period $=$'+str("{:.2f}".format(round(t_text[0],2))),fontsize=6)
                
            self.ax3.legend(fontsize=6)
            return self.v_x_scat,self.v1_x_scat,self.v2_x_scat,self.v_y_scat,self.v1_y_scat,self.v2_y_scat,self.v_z_scat,self.v1_z_scat,self.v2_z_scat,self.v_obs_x_scat,self.v1_obs_x_scat,self.v2_obs_x_scat,self.v_obs_y_scat,self.v1_obs_y_scat,self.v2_obs_y_scat,self.v_obs_z_scat,self.v1_obs_z_scat,self.v2_obs_z_scat,self.t_text
    
        def data_stream(self):
            
            i=0
            while True:
                v_x=sample_orbit.v_xyz.T[0][i].value/np.sqrt(-2.*sample_orbit.energy[i].value)
                v1_x=sample_orbit.v1_xyz.T[0][i].value/np.sqrt(-2.*sample_orbit.energy[i].value)
                v2_x=sample_orbit.v2_xyz.T[0][i].value/np.sqrt(-2.*sample_orbit.energy[i].value)
                v_y=sample_orbit.v_xyz.T[1][i].value/np.sqrt(-2.*sample_orbit.energy[i].value)
                v1_y=sample_orbit.v1_xyz.T[1][i].value/np.sqrt(-2.*sample_orbit.energy[i].value)
                v2_y=sample_orbit.v2_xyz.T[1][i].value/np.sqrt(-2.*sample_orbit.energy[i].value)
                v_z=sample_orbit.v_xyz.T[2][i].value/np.sqrt(-2.*sample_orbit.energy[i].value)
                v1_z=sample_orbit.v1_xyz.T[2][i].value/np.sqrt(-2.*sample_orbit.energy[i].value)
                v2_z=sample_orbit.v2_xyz.T[2][i].value/np.sqrt(-2.*sample_orbit.energy[i].value)

                v_obs_x=sample_orbit.v_obs_xyz.T[0].value[i]/np.sqrt(-2.*sample_orbit.energy[i].value)
                v1_obs_x=sample_orbit.v1_obs_xyz.T[0].value[i]/np.sqrt(-2.*sample_orbit.energy[i].value)
                v2_obs_x=sample_orbit.v2_obs_xyz.T[0].value[i]/np.sqrt(-2.*sample_orbit.energy[i].value)
                v_obs_y=sample_orbit.v_obs_xyz.T[1].value[i]/np.sqrt(-2.*sample_orbit.energy[i].value)
                v1_obs_y=sample_orbit.v1_obs_xyz.T[1].value[i]/np.sqrt(-2.*sample_orbit.energy[i].value)
                v2_obs_y=sample_orbit.v2_obs_xyz.T[1].value[i]/np.sqrt(-2.*sample_orbit.energy[i].value)
                v_obs_z=sample_orbit.v_obs_xyz.T[2].value[i]/np.sqrt(-2.*sample_orbit.energy[i].value)
                v1_obs_z=sample_orbit.v1_obs_xyz.T[2].value[i]/np.sqrt(-2.*sample_orbit.energy[i].value)
                v2_obs_z=sample_orbit.v2_obs_xyz.T[2].value[i]/np.sqrt(-2.*sample_orbit.energy[i].value)
            
                t=sample_orbit.time[i].value/sample_orbit.period[i].value
                i+=1
                if i==len(sample_orbit.r_xyz):#infinite loop
                    i=0
                yield np.c_[t,v_x,v1_x,v2_x,v_y,v1_y,v2_y,v_z,v1_z,v2_z,v_obs_x,v1_obs_x,v2_obs_x,v_obs_y,v1_obs_y,v2_obs_y,v_obs_z,v1_obs_z,v2_obs_z,t]

        def update(self, i):
            
            data = next(self.stream)

            self.v_x_scat.set_offsets([[data[0][0],data[0][1]]])
            self.v1_x_scat.set_offsets([[data[0][0],data[0][2]]])
            self.v2_x_scat.set_offsets([[data[0][0],data[0][3]]])
            self.v_y_scat.set_offsets([[data[0][0],data[0][4]]])
            self.v1_y_scat.set_offsets([[data[0][0],data[0][5]]])
            self.v2_y_scat.set_offsets([[data[0][0],data[0][6]]])
            self.v_z_scat.set_offsets([[data[0][0],data[0][7]]])
            self.v1_z_scat.set_offsets([[data[0][0],data[0][8]]])
            self.v2_z_scat.set_offsets([[data[0][0],data[0][9]]])
            
            self.v_obs_x_scat.set_offsets([[data[0][0],data[0][10]]])
            self.v1_obs_x_scat.set_offsets([[data[0][0],data[0][11]]])
            self.v2_obs_x_scat.set_offsets([[data[0][0],data[0][12]]])
            self.v_obs_y_scat.set_offsets([[data[0][0],data[0][13]]])
            self.v1_obs_y_scat.set_offsets([[data[0][0],data[0][14]]])
            self.v2_obs_y_scat.set_offsets([[data[0][0],data[0][15]]])
            self.v_obs_z_scat.set_offsets([[data[0][0],data[0][16]]])
            self.v1_obs_z_scat.set_offsets([[data[0][0],data[0][17]]])
            self.v2_obs_z_scat.set_offsets([[data[0][0],data[0][18]]])
        
            self.t_text.set_text('time / period='+str("{:.2f}".format(round(data[0][0],2))))

            #return the updated artist
            return self.v_x_scat,self.v1_x_scat,self.v2_x_scat,self.v_y_scat,self.v1_y_scat,self.v2_y_scat,self.v_z_scat,self.v1_z_scat,self.v2_z_scat,self.v_obs_x_scat,self.v1_obs_x_scat,self.v2_obs_x_scat,self.v_obs_y_scat,self.v1_obs_y_scat,self.v2_obs_y_scat,self.v_obs_z_scat,self.v1_obs_z_scat,self.v2_obs_z_scat,self.t_text
    a=animate_v()
    plt.show()

def add_binaries(object_xyz,mass_primary,**params):#mass is mass_primary+mass_secondary
    
    class r2d_with_binaries:    
        def __init__(self,r_xyz=None,mass=None,item=None,companion=None,binary_model=None):
            self.r_xyz=r_xyz #AU
            self.mass=mass
            self.item=item
            self.companion=companion
            self.binary_model=binary_model

    if not 'binary_model' in params:
        params['binary_model']='user'
    if not 'f_binary' in params:
        params['f_binary']=1.

    n_object=len(object_xyz)
    
    is_binary=np.zeros(n_object,dtype='bool')
    is_binary[np.random.uniform(size=n_object,low=0.,high=1.)<=params['f_binary']]=True
    
    n_binary=is_binary.sum()
    n_single=n_object-n_binary

    m_min=0.1 #Msun, minimum mass of star in binary system
    
    if params['binary_model']=='Raghavan2010':

        q=np.random.uniform(size=n_object,low=m_min/mass_primary.value,high=1.) #array of m_secondary / m_primary, sampled from uniform distribution subject to constraint M_secondary > M_min
        period=10.**sample_normal_truncated(size=n_object,loc=5.03,scale=2.28,min_value=-np.inf,max_value=np.inf)/364.25*u.yr #array of orbital period (years), sampled from truncated log-normal distribution
        eccentricity=np.random.uniform(size=n_object,low=0.,high=1.)
        #eccentricity=10.**sample_normal_truncated(size=n_binary,loc=-0.3,scale=1.,min_value=-np.inf,max_value=0.) #array of orbital eccentricity, sampled from truncated log-normal distribution
        eccentricity[period*365.24<12.*u.day]=0. #eccentricity=0 for P<12 days

    elif params['binary_model']=='DM91':
        q=sample_normal_truncated(size=n_object,loc=0.23,scale=0.42,min_value=m_min/mass_primary.value,max_value=1.)
        period=10.**sample_normal_truncated(size=n_object,loc=4.8,scale=2.3,min_value=-np.inf,max_value=np.inf)/364.25*u.yr #array of orbital period (years), sampled from truncated log-normal distribution
        eccentricity=sample_normal_truncated(size=n_object,loc=0.31,scale=0.17,min_value=0.,max_value=1.)
        long_period=np.where(period*365.24>1000.*u.day)[0]
        eccentricity_thermal=sample_thermal(size=len(long_period)) #sample thermal distribution for long periods
        eccentricity[long_period]=eccentricity_thermal
        eccentricity[period*365.24<12.*u.day]=0. #eccentricity=0 for P<12 days

    else:
        q=params['q']
        period=params['period']
        eccentricity=params['eccentricity']
        
    f_period=np.random.uniform(size=n_object,low=0.,high=1.) #array of orbital phase, time / period
    inclination=sample_inclination(size=n_object)*u.rad #array of inclination angle (radians), inclination=0 for observer along +z axis, inclination=pi/2 for observer in xy plane, allowed from 0 to 2*pi to allow for full range of parity.
    longitude=np.random.uniform(size=n_object,low=0,high=2.*np.pi)*u.rad #array of longitude of ascending node (radians), longitude=0 if observer is along +x axis, longitude=pi/2 if observer is along +y axis

    orbit_snapshot=sample_orbit_2body(f_period,period=period,eccentricity=eccentricity,mass_primary=mass_primary,mass_ratio=q,longitude=longitude,inclination=inclination)

    r_xyz=np.zeros((n_object+n_binary,3))*object_xyz[0].unit
    mass=np.zeros(n_object+n_binary)*orbit_snapshot.mass_primary[0].unit
    item=np.zeros(n_object+n_binary,dtype='int')
    companion=np.zeros(n_object+n_binary,dtype='int')

    j=0
    
    for i in range(0,len(object_xyz)):
        
        if is_binary[i]:

            if np.abs(1.-mass_primary[i]/orbit_snapshot.mass_primary[i])>0.01:
                raise ValueError('problem with binary mass samples!')
            r1=(object_xyz[i]+orbit_snapshot.r1_obs_xyz[i]).to(object_xyz[0].unit)
            r2=(object_xyz[i]+orbit_snapshot.r2_obs_xyz[i]).to(object_xyz[0].unit)
            
            r_xyz.value[j]=r1.value
            mass.value[j]=orbit_snapshot.mass_primary.value[i]
            item[j]=1
            companion[j]=j+1
            j+=1
            
            r_xyz.value[j]=r2.value
            mass.value[j]=orbit_snapshot.mass_secondary.value[i]
            item[j]=2
            companion[j]=j-1
            j+=1
            
        else:
            
            r_xyz.value[j]=object_xyz[i].value
            mass.value[j]=mass_primary.value[i]
            item[j]=0
            companion[j]=-999
            j+=1
            
    return r2d_with_binaries(r_xyz=np.array(r_xyz)*r1.unit,mass=np.array(mass),item=np.array(item),companion=np.array(companion,dtype=int),binary_model=params['binary_model'])


