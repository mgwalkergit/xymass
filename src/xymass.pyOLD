import numpy as np
import scipy
import scipy.optimize
import scipy.special
import warnings

def sample_r2d(size,model,**params):#samples from flattened plummer, exponential, or (not flattened) uniform 2d distributions
    
    class r2d:
        def __init__(self,r_ell=None,x=None,y=None,ellipticity=None,position_angle=None,r_scale=None,model=None,alpha=None,beta=None,gamma=None,func=None):
            self.r_ell=r_ell
            self.x=x
            self.y=y
            self.ellipticity=ellipticity
            self.position_angle=position_angle
            self.r_scale=r_scale
            self.model=model
            self.alpha=alpha
            self.beta=beta
            self.gamma=gamma
            self.func=func

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
        def func(x):
            return bigsigma0/(1+x**2)**2
        r=np.sqrt(uni/(1.-uni))#elliptical radius
        return r2d(r_ell=r*params['r_scale'],x=r*flat_x*params['r_scale'],y=r*flat_y*params['r_scale'],ellipticity=params['ellipticity'],position_angle=params['position_angle'],r_scale=params['r_scale'],model=model,func=func)

    if model=='exp':
        bigsigma0=size/2/np.pi/params['r_scale']**2
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
        return r2d(r_ell=r*params['r_scale'],x=r*flat_x*params['r_scale'],y=r*flat_y*params['r_scale'],ellipticity=params['ellipticity'],position_angle=params['position_angle'],r_scale=params['r_scale'],model=model,func=func)

    if model=='2bg':
        bigsigma0=size*(params['beta']-3)*scipy.special.gamma((params['beta']-params['gamma'])/2)/4/np.sqrt(np.pi)/scipy.special.gamma((3-params['gamma'])/2)/scipy.special.gamma(params['beta']/2)/params['r_scale']**2
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
        return r2d(r_ell=r*params['r_scale'],x=r*flat_x*params['r_scale'],y=r*flat_y*params['r_scale'],ellipticity=params['ellipticity'],position_angle=params['position_angle'],r_scale=params['r_scale'],model=model,beta=params['beta'],gamma=params['gamma'],func=func)
    
    if model=='uni':
        bigsigma0=size/np.pi/params['r_scale']**2
        def func(x):
            return bigsigma0*x/x
        r=np.sqrt(uni)#elliptical radius (can in practice be elliptical if nonzero ellipticity is specified)
        return r2d(r_ell=r*params['r_scale'],x=r*flat_x*params['r_scale'],y=r*flat_y*params['r_scale'],ellipticity=params['ellipticity'],position_angle=params['position_angle'],r_scale=params['r_scale'],model=model,func=func)


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
        params['m_min']=0.03
    if not 'm_max' in params:
        params['m_max']=150.
        
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
            raise ValueError ('problem with parameters in Kroupa IMF')
        if params['m1_break']<params['m_min']:
            raise ValueError ('problem with parameters in Kroupa IMF')
        if params['m2_break']>params['m_max']:
            raise ValueError ('problem with parameters in Kroupa IMF')        
            
        mass=[]
        piece1=(params['m1_break']**(1.-params['alpha1'])-params['m_min']**(1.-params['alpha1']))/(1.-params['alpha1'])
        piece2=(params['m2_break']**(1.-params['alpha2'])-params['m1_break']**(1.-params['alpha2']))/(1.-params['alpha2'])
        piece3=(params['m_max']**(1.-params['alpha3'])-params['m2_break']**(1.-params['alpha3']))/(1.-params['alpha3'])
        #get normalization constant for each of three pieces
        k2_over_k1=params['m1_break']**(params['alpha2']-params['alpha1'])
        k3_over_k2=params['m2_break']**(params['alpha3']-params['alpha2'])
        
        k1=1./(piece1+piece2*k2_over_k1+piece3*k3_over_k2*k2_over_k1)#sample size normalized to 1
        k2=k1*k2_over_k1
        k3=k2*k3_over_k2

        #get fraction of sample within each piece
        f1=k1*piece1
        f2=k2*piece2
        f3=k3*piece3

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
        
        for i in range(0,len(ran1)):
            if ran1[i]<f1:
                mass.append((params['m_min']**(1.-params['alpha1'])+ran1[i]*(1.-params['alpha1'])/k1)**(1./(1.-params['alpha1'])))
            elif ((ran1[i]>=f1)&(ran1[i]<(f1+f2))):
                mass.append((params['m1_break']**(1.-params['alpha2'])+(1.-params['alpha2'])/k2*(ran1[i]-f1))**(1./(1.-params['alpha2'])))
            elif ran1[i]>=(f1+f2):
                mass.append((params['m2_break']**(1.-params['alpha3'])+(1.-params['alpha3'])/k3*(ran1[i]-f1-f2))**(1./(1.-params['alpha3'])))
            else:
                raise ValueError('something wrong in sampling Kroupa IMF')
        mass=np.array(mass)
        return imf(model=model,mass=mass,alpha1=params['alpha1'],alpha2=params['alpha2'],alpha3=params['alpha3'],m1_break=params['m1_break'],m2_break=params['m2_break'],m_min=params['m_min'],m_max=params['m_max'],k1=k1,k2=k2,k3=k3,func=kroupa_func)

    if model=='bpl':#sample from BPL IMF, 2 separate power laws with indices -alpha1, -alpha2, break mass at m_break

        if not 'alpha1' in params:
            params['alpha1']=1.3
        if not 'alpha2' in params:
            params['alpha2']=2.3
        if not 'm_break' in params:
            params['m_break']=0.5
        
        if params['m_break']<params['m_min']:
            raise ValueError ('problem with parameters in BPL IMF')
        if params['m_break']>params['m_max']:
            raise ValueError ('problem with parameters in BPL IMF')
        
        mass=[]
        piece1=(params['m_break']**(1.-params['alpha1'])-params['m_min']**(1.-params['alpha1']))/(1.-params['alpha1'])
        piece2=(params['m_max']**(1.-params['alpha2'])-params['m_break']**(1.-params['alpha2']))/(1.-params['alpha2'])
        #get normalization constant for each of three pieces
        k2_over_k1=params['m_break']**(params['alpha2']-params['alpha1'])
        
        k1=1./(piece1+piece2*k2_over_k1)#sample size normalized to 1
        k2=k1*k2_over_k1

        #get fraction of sample within each piece
        f1=k1*piece1
        f2=k2*piece2

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
                
            elif ((type(x) is float)|(type(x) is int)|(type(x) is np.float64)):
                if x<params['m_break']:
                    val=k1*x**-params['alpha1']
                elif x>=params['m_break']:
                    val=k2*x**-params['alpha2']
                else:
                    raise ValueError('problem in bpl_func')
            else:
                raise ValueError('type error in bpl_func')
                
            return val
        
        for i in range(0,len(ran1)):
            if ran1[i]<f1:
                mass.append((params['m_min']**(1.-params['alpha1'])+ran1[i]*(1.-params['alpha1'])/k1)**(1./(1.-params['alpha1'])))
            elif ran1[i]>=f1:
                mass.append((params['m_break']**(1.-params['alpha2'])+(1.-params['alpha2'])/k2*(ran1[i]-f1))**(1./(1.-params['alpha2'])))
            else:
                raise ValueError('something wrong in sampling BPL IMF')
        mass=np.array(mass)
        return imf(model=model,mass=mass,alpha1=params['alpha1'],alpha2=params['alpha2'],m_break=params['m_break'],m_min=params['m_min'],m_max=params['m_max'],k1=k1,k2=k2,func=bpl_func)

