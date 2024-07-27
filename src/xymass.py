import sampler
import numpy as np
import scipy
import scipy.optimize
import scipy.special
import warnings
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
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

        r=sampler.plum(size)#elliptical radius/r_scale

        return r2d(r_ell=r*params['r_scale'],x=r*flat_x*params['r_scale'],y=r*flat_y*params['r_scale'],r_xyz=np.c_[r*flat_x*params['r_scale'],r*flat_y*params['r_scale'],np.zeros(len(r),dtype=float)],ellipticity=params['ellipticity'],position_angle=params['position_angle'],r_scale=params['r_scale'],model=model,func=func,rhalf_2d=rhalf_2d)

    if model=='exp':
        
        bigsigma0=size/2/np.pi/params['r_scale']**2
        rhalf_2d=1.67835*params['r_scale']

        def func(x):
            return bigsigma0*np.exp(-x)

        r=sampler.exp(size) #elliptical radius / r_scale

        return r2d(r_ell=r*params['r_scale'],x=r*flat_x*params['r_scale'],y=r*flat_y*params['r_scale'],r_xyz=np.c_[r*flat_x*params['r_scale'],r*flat_y*params['r_scale'],np.zeros(len(r),dtype=float)],ellipticity=params['ellipticity'],position_angle=params['position_angle'],r_scale=params['r_scale'],model=model,func=func,rhalf_2d=rhalf_2d)

    if model=='a2bg':
        
        bigsigma0=size*(params['beta']-3)*scipy.special.gamma((params['beta']-params['gamma'])/2)/4/np.sqrt(np.pi)/scipy.special.gamma((3-params['gamma'])/2)/scipy.special.gamma(params['beta']/2)/params['r_scale']**2

        def rootfind_a2bg_2d(x,beta,gamma):
            return 0.5-np.sqrt(np.pi)*scipy.special.gamma((beta-gamma)/2)/2/scipy.special.gamma(beta/2)/scipy.special.gamma((3-gamma)/2)*x**(3-beta)*scipy.special.hyp2f1((beta-3)/2,(beta-gamma)/2,beta/2,-1/x**2)

        low0=1.e-10
        high0=1.e+10
        rhalf_2d=params['r_scale']*scipy.optimize.brentq(rootfind_a2bg_2d,low0,high0,args=(params['beta'],params['gamma']),xtol=1.e-12,rtol=1.e-6,maxiter=1000,full_output=False,disp=True)

        def func(x):
            return bigsigma0*x**(1-params['beta'])*scipy.special.hyp2f1((params['beta']-1)/2,(params['beta']-params['gamma'])/2,params['beta']/2,-1/x**2)            

        r=sampler.a2bg(size,params['beta'],params['gamma'])#elliptical radius / r_scale

        return r2d(r_ell=r*params['r_scale'],x=r*flat_x*params['r_scale'],y=r*flat_y*params['r_scale'],r_xyz=np.c_[r*flat_x*params['r_scale'],r*flat_y*params['r_scale'],np.zeros(len(r),dtype=float)],ellipticity=params['ellipticity'],position_angle=params['position_angle'],r_scale=params['r_scale'],model=model,beta=params['beta'],gamma=params['gamma'],func=func,rhalf_2d=rhalf_2d)
    
    if model=='uni':
        bigsigma0=size/np.pi/params['r_scale']**2
        rhalf_2d=np.sqrt(0.5)*params['r_scale']
        
        def func(x):
            return bigsigma0*x/x
        
        r=sampler.uni(size)#elliptical radius (can in practice be elliptical if nonzero ellipticity is specified) / r_scale

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
        
    if model=='salpeter':

        if not 'alpha' in params:
            params['alpha']=2.3
        
        k_salpeter=(1.-params['alpha'])/(params['m_max']**(1.-params['alpha'])-params['m_min']**(1.-params['alpha']))
        def salpeter_func(x):
            return k_salpeter*x**-params['alpha']
            
        mass=sampler.pl(size,params['m_min'],params['m_max'],params['alpha'])
        
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
            
        mass=sampler.lognormal(size,params['m_min'],params['m_max'],params['mean'],params['std'])
        
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
        
        #get normalization constant for each of three pieces
        k2_over_k1=params['m1_break']**(params['alpha2']-params['alpha1'])
        k3_over_k2=params['m2_break']**(params['alpha3']-params['alpha2'])
        
        mass,k1,k2,k3=sampler.kroupa(size,params['m_min'],params['m_max'],params['alpha1'],params['alpha2'],params['alpha3'],params['m1_break'],params['m2_break'])
                            
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
            
        #get normalization constant for each of three pieces
        k2_over_k1=params['m_break']**(params['alpha2']-params['alpha1'])
        
        mass,k1,k2=sampler.bpl(size,params['m_min'],params['m_max'],params['alpha1'],params['alpha2'],params['m_break'])
        
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

    return sampler.normal_truncated(params['size'],params['min_value'],params['max_value'],params['loc'],params['scale'])

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

def add_binaries_physical(object_xyz,mass_primary,**params):
    
    class r2d_with_binaries:    
        def __init__(self,r_xyz=None,mass=None,item=None,companion=None,binary_model=None):
            self.r_xyz=r_xyz
            self.mass=mass
            self.item=item
            self.companion=companion
            self.binary_model=binary_model

    if not 'binary_model' in params:
        params['binary_model']='user'
    if not 'f_binary' in params:
        params['f_binary']=1.
    if not 'm_min' in params:
        params['m_min']=np.min(mass_primary.value)
    if (('mass_secondary' in params)&('mass_ratio' in params)):
        raise ValueError('cannot specify both mass_seconary and mass_ratio')
    if 'mass_secondary' in params:
        params['mass_ratio']=params['mass_secondary']/mass_primary
    elif 'mass_ratio' in params:
        params['mass_secondary']=mass_primary*params['mass_ratio']
        
    n_object=len(object_xyz)
    
    is_binary=np.zeros(n_object,dtype='bool')
    is_binary[np.random.uniform(size=n_object,low=0.,high=1.)<=params['f_binary']]=True
    
    n_binary=is_binary.sum()
    n_single=n_object-n_binary

    if params['binary_model']=='Raghavan2010':

        mass_ratio=np.random.uniform(size=n_object,low=params['m_min']/mass_primary.value,high=1.) #array of m_secondary / m_primary, sampled from uniform distribution subject to constraint M_secondary > M_min
        period=10.**sample_normal_truncated(size=n_object,loc=5.03,scale=2.28,min_value=-np.inf,max_value=np.inf)/364.25*u.yr #array of orbital period (years), sampled from truncated log-normal distribution
        eccentricity=np.random.uniform(size=n_object,low=0.,high=1.)
        #eccentricity=10.**sample_normal_truncated(size=n_binary,loc=-0.3,scale=1.,min_value=-np.inf,max_value=0.) #array of orbital eccentricity, sampled from truncated log-normal distribution
        eccentricity[period*365.24<12.*u.day]=0. #eccentricity=0 for P<12 days

    elif params['binary_model']=='DM91':
        mass_ratio=sample_normal_truncated(size=n_object,loc=0.23,scale=0.42,min_value=params['m_min']/mass_primary.value,max_value=1.)
        period=10.**sample_normal_truncated(size=n_object,loc=4.8,scale=2.3,min_value=-np.inf,max_value=np.inf)/364.25*u.yr #array of orbital period (years), sampled from truncated log-normal distribution
        eccentricity=sample_normal_truncated(size=n_object,loc=0.31,scale=0.17,min_value=0.,max_value=1.)
        long_period=np.where(period*365.24>1000.*u.day)[0]
        eccentricity_thermal=sampler.uni(size=len(long_period)) #sample thermal distribution for long periods, this is equivalent to sampling radial coordinate of uniform 2D distribution
        eccentricity[long_period]=eccentricity_thermal
        eccentricity[period*365.24<12.*u.day]=0. #eccentricity=0 for P<12 days

    else:
        mass_ratio=params['mass_ratio']
        period=params['period']
        eccentricity=params['eccentricity']
        
    f_period=np.random.uniform(size=n_object,low=0.,high=1.) #array of orbital phase, time / period
    inclination=sample_inclination(size=n_object)*u.rad #array of inclination angle (radians), inclination=0 for observer along +z axis, inclination=pi/2 for observer in xy plane, allowed from 0 to 2*pi to allow for full range of parity.
    longitude=np.random.uniform(size=n_object,low=0,high=2.*np.pi)*u.rad #array of longitude of ascending node (radians), longitude=0 if observer is along +x axis, longitude=pi/2 if observer is along +y axis

    orbit_snapshot=sample_orbit_2body(f_period,period=period,eccentricity=eccentricity,mass_primary=mass_primary,mass_ratio=mass_ratio,longitude=longitude,inclination=inclination)

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

def add_binaries_func(object_xyz,**params):
    
    class r2d_with_binaries_func:    
        def __init__(self,r_xyz=None,mass=None,item=None,companion=None,separation_func=None,projected=None):
            self.r_xyz=r_xyz
            self.mass=mass
            self.item=item
            self.companion=companion
            self.separation_func=separation_func
            self.projected=projected

    if not 'f_binary' in params:
        params['f_binary']=1.
    if (('mass_secondary' in params)&('mass_ratio' in params)):
        raise ValueError('cannot specify both mass_seconary and mass_ratio')
    if 'mass_secondary' in params:
        params['mass_ratio']=params['mass_secondary']/params['mass_primary']
    elif 'mass_ratio' in params:
        params['mass_secondary']=params['mass_primary']*params['mass_ratio']
        
    n_object=len(object_xyz)
    
    is_binary=np.zeros(n_object,dtype='bool')
    is_binary[np.random.uniform(size=n_object,low=0.,high=1.)<=params['f_binary']]=True
    
    n_binary=is_binary.sum()
    n_single=n_object-n_binary

    if params['separation_func']=='pl':
        r=sampler.pl(len(object_xyz),params['s_min'].to(u.AU).value,params['s_max'].to(u.AU).value,params['alpha'])*params['s_min'].unit
        
    if params['separation_func']=='bpl':
        r=sampler.bpl(len(object_xyz),params['s_min'].to(u.AU).value,params['s_max'].to(u.AU).value,params['alpha1'],params['alpha2'],params['s_break'].to(u.AU).value)[0]*params['s_min'].unit

    if params['separation_func']=='lognormal':
        r=10.**sampler.normal_truncated(len(object_xyz),np.log10(params['s_min'].to(u.AU).value),np.log10(params['s_max'].to(u.AU).value),params['loc'],params['scale'])*params['s_min'].unit
                                        
    longitude=np.random.uniform(size=n_object,low=0,high=2.*np.pi)*u.rad
    if params['projected']:
        inclination=np.zeros(len(object_xyz),dtype=float)*u.rad #if separation function is projected, view binary orbit face-on
    else:
        inclination=sample_inclination(size=len(object_xyz))*u.rad
        
    x=r#*np.cos(theta) #x component of separation vector in orbital plane, effectively assume theta=0 (pericenter, if this were orbit calculation)
    y=r-r#sep*np.sin(theta) # y component of separation vector in orbital plane, effectively assume theta=0 (pericenter, if this were orbit calculation)
    z=r-r
        
    trans1,trans2=-params['mass_ratio']/(1.+params['mass_ratio']),1./(1.+params['mass_ratio'])    
    x1,y1=trans1*x,trans1*y
    x2,y2=trans2*x,trans2*y

    r_xyz=np.array((x,y,z)).T*r.unit
    r1_xyz=np.array((x1,y1,z)).T*r.unit
    r2_xyz=np.array((x2,y2,z)).T*r.unit

    r_obs_xyz=np.zeros(np.shape(r_xyz))*r.unit
    r1_obs_xyz=np.zeros(np.shape(r_xyz))*r.unit
    r2_obs_xyz=np.zeros(np.shape(r_xyz))*r.unit
    
    rot_matrix=[]
        
    rot_alpha=longitude.to(u.rad).value #rotation about z axis, in direction of arc from +x to +y (radians)
    rot_beta=0. #rotation about y axis, in direction of arc from +z to +x (radians)
    rot_gamma=inclination.to(u.rad).value #rotation about x axis, in direction of arc from +y to +z (radians)

    for i in range(0,len(longitude)):
        rot_matrix.append(get_rot_matrix(rot_alpha[i],rot_beta,rot_gamma[i]))

    for i in range(0,len(rot_matrix)):
        r_obs_xyz.value[i]=rot_matrix[i].apply(r_xyz.value[i])
        r1_obs_xyz.value[i]=rot_matrix[i].apply(r1_xyz.value[i])
        r2_obs_xyz.value[i]=rot_matrix[i].apply(r2_xyz.value[i])    

    r_xyz=np.zeros((n_object+n_binary,3))*object_xyz[0].unit
    mass=np.zeros(n_object+n_binary)*params['mass_primary'][0].unit
    item=np.zeros(n_object+n_binary,dtype='int')
    companion=np.zeros(n_object+n_binary,dtype='int')

    j=0
    
    for i in range(0,len(object_xyz)):
        
        if is_binary[i]:

            if np.abs(1.-params['mass_primary'][i]/params['mass_primary'][i])>0.01:
                raise ValueError('problem with binary mass samples!')
            r1=(object_xyz[i]+r1_obs_xyz[i]).to(object_xyz[0].unit)
            r2=(object_xyz[i]+r2_obs_xyz[i]).to(object_xyz[0].unit)
            
            r_xyz.value[j]=r1.value
            mass.value[j]=params['mass_primary'].value[i]
            item[j]=1
            companion[j]=j+1
            j+=1
            
            r_xyz.value[j]=r2.value
            mass.value[j]=params['mass_secondary'].value[i]
            item[j]=2
            companion[j]=j-1
            j+=1
            
        else:
            
            r_xyz.value[j]=object_xyz[i].value
            mass.value[j]=params['mass_primary'].value[i]
            item[j]=0
            companion[j]=-999
            j+=1
            
    return r2d_with_binaries_func(r_xyz=np.array(r_xyz)*r1.unit,mass=np.array(mass),item=np.array(item),companion=np.array(companion,dtype=int),separation_func=params['separation_func'],projected=params['projected'])

def binary_blend(r2d_wb,mag,Mbol_sun,resolution_limit_physical,**params):
    
    class r2d_binary_blend:
        def __init__(self,r_xy=None,mag=None,blend=None):
            self.r_xy=r_xy
            self.mag=mag
            self.blend=blend

    xy=np.c_[r2d_wb.r_xyz.T[0],r2d_wb.r_xyz.T[1]]
    mass=np.c_[r2d_wb.mass,r2d_wb.mass]

    tree=scipy.spatial.KDTree(xy)
    #tree_single=scipy.spatial.KDTree(xy[r2d_wb.item==0])
    #tree_binary=scipy.spatial.KDTree(xy[r2d_wb.item>0])
    query=tree.query(xy,k=2)
    #query_single=tree_single.query(xy[r2d_wb.item==0],k=2)
    #query_binary=tree_binary.query(xy[r2d_wb.item>0],k=2)
    nn=query[0].T[1]
    #nn_single=query_single[0].T[1]
    #nn_binary=query_binary[0].T[1]

    nn_partner=query[1].T[1]

    lum=10.**(-(mag-Mbol_sun)/2.5)
    
    unresolved=np.where(nn*xy.unit<resolution_limit_physical)[0]
    
    xy[unresolved]=(mass[unresolved]*xy[unresolved]+mass[nn_partner[unresolved]]*xy[nn_partner[unresolved]])/(mass[unresolved]+mass[nn_partner[unresolved]]) #replace position with mass-weighted mean position of unresolved partners
    lum[unresolved]=lum[unresolved]+lum[nn_partner[unresolved]] #replace luminosity with sum of luminosities of unresolved partners
    mag[unresolved]=Mbol_sun-2.5*np.log10(lum[unresolved]) #magnitude corresponding to sum of luminoisities

    keep=np.full(len(xy),True,dtype='bool')
    blend=np.full(len(xy),False,dtype='bool')
    blend[unresolved]=True
    for i in range(0,len(unresolved)):
        if keep[unresolved[i]]:
            keep[nn_partner[unresolved[i]]]=False
            
    return r2d_binary_blend(r_xy=xy[keep],mag=mag[keep],blend=blend[keep])

