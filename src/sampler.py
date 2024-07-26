import numpy as np
import scipy
import scipy.optimize
import scipy.special

low0=1.e-10
high0=1.e10

def plum(size):#sample radial coordinate from plummer profile
    ran=np.random.uniform(low=0.,high=1.,size=size)
    return np.sqrt(ran/(1.-ran))

def exp(size):#sample radial coordinate from exponential profile
    ran=np.random.uniform(low=0.,high=1.,size=size)
    def findx_exp(x,ran):
        return 1.-(1.+x)*np.exp(-x)-ran
    x=np.zeros(size,dtype='float')
    for i in range(0,len(ran)):
        x[i]=scipy.optimize.brentq(findx_exp,low0,high0,args=ran[i],xtol=1.e-12,rtol=1.e-6,maxiter=100,full_output=False,disp=True)
    return x

def a2bg(size,beta,gamma):#sample radial coordinate from alpha/beta/gamma profile with alpha=2
    ran=np.random.uniform(low=0.,high=1.,size=size)
    def findx_a2bg(arg,rand,beta,gamma):
        return 1-np.sqrt(np.pi)/2*scipy.special.gamma((beta-gamma)/2)/scipy.special.gamma(beta/2)/scipy.special.gamma((3-gamma)/2)*arg**(3-beta)*scipy.special.hyp2f1((beta-3)/2,(beta-gamma)/2,beta/2,-1/arg**2)-rand
    
    x=np.zeros(size,dtype='float')
    for i in range(0,len(ran)):
        x[i]=scipy.optimize.brentq(findx_a2bg,low0,high0,args=(ran[i],beta,gamma),xtol=1.e-12,rtol=1.e-6,maxiter=100,full_output=False,disp=True)
    return x

def uni(size):#sample radial coordinate from uniform (2D) distribution
    ran=np.random.uniform(low=0.,high=1.,size=size)
    return np.sqrt(ran)

def pl(size,x_min,x_max,alpha):#sample coordinate from power law
    ran=np.random.uniform(low=0.,high=1.,size=size)
    return (x_min**(1.-alpha)+ran*(x_max**(1.-alpha)-x_min**(1.-alpha)))**(1./(1.-alpha))

def lognormal(size,x_min,x_max,loc,scale):
    ran=np.random.uniform(low=0.,high=1.,size=size)
    ntotnorm=scipy.special.erf((np.log10(loc)*np.log(10.)-np.log(x_min))/np.sqrt(2.)/np.log(10.)/scale)-scipy.special.erf((np.log10(loc)*np.log(10.)-np.log(x_max))/np.sqrt(2.)/np.log(10.)/scale)
    erf=scipy.special.erf((np.log10(loc)*np.log(10.)-np.log(x_min))/np.sqrt(2.)/np.log(10.)/scale)-ran*ntotnorm
    return np.exp(np.log10(loc)*np.log(10.)-np.sqrt(2.)*np.log(10.)*scale*scipy.special.erfinv(erf))

def kroupa(size,x_min,x_max,alpha1,alpha2,alpha3,x1,x2):
    ran=np.random.uniform(low=0.,high=1.,size=size)
    x=np.zeros(size,dtype=float)

    #get normalization constant for each of three pieces
    k2_over_k1=x1**(alpha2-alpha1)
    k3_over_k2=x2**(alpha3-alpha2)
        
    if x_min<x1:
            
        if x_max>x2:
                
            piece1=(x1**(1.-alpha1)-x_min**(1.-alpha1))/(1.-alpha1)
            piece2=(x2**(1.-alpha2)-x1**(1.-alpha2))/(1.-alpha2)
            piece3=(x_max**(1.-alpha3)-x2**(1.-alpha3))/(1.-alpha3)

            x0_2=x1
            x0_3=x2
                                
        if ((x1<=x_max)&(x_max<=x2)):
                
            piece1=(x1**(1.-alpha1)-x_min**(1.-alpha1))/(1.-alpha1)
            piece2=(x_max**(1.-alpha2)-x1**(1.-alpha2))/(1.-alpha2)
            piece3=0.

            x0_2=x1
            x0_3=x2
                                
        if x1>x_max:
                
            piece1=(x_max**(1.-alpha1)-x_min**(1.-alpha1))/(1.-alpha1)
            piece2=0.
            piece3=0.
            
            x0_2=x1
            x0_3=x2
                
    if ((x1<=x_min)&(x_min<=x2)):
            
        if x_max>x2:
                
            piece1=0.
            piece2=(x2**(1.-alpha2)-x_min**(1.-alpha2))/(1.-alpha2)
            piece3=(x_max**(1.-alpha3)-x2**(1.-alpha3))/(1.-alpha3)
                
            x0_2=x_min
            x0_3=x2
            
        if ((x1<=x_max)&(x_max<=x2)):
                
            piece1=0.
            piece2=(x_max**(1.-alpha2)-x_min**(1.-alpha2))/(1.-alpha2)
            piece3=0.
                
            x0_2=x_min
            x0_3=x2
                
    if x_min>x2:
            
        if x_max>x2:
                
            piece1=0.
            piece2=0.
            piece3=(x_max**(1.-alpha3)-x_min**(1.-alpha3))/(1.-alpha3)
                
            x0_2=x1
            x0_3=x_min
                
    k1=1./(piece1+piece2*k2_over_k1+piece3*k3_over_k2*k2_over_k1)#sample size normalized to 1
    k2=k1*k2_over_k1
    k3=k2*k3_over_k2

    #get fraction of sample within each piece
    f1=k1*piece1
    f2=k2*piece2
    f3=k3*piece3

    first=np.where(ran<f1)[0]
    second=np.where((ran>=f1)&(ran<f1+f2))[0]
    third=np.where(ran>=f1+f2)[0]
    bad=np.where(ran>f1+f2+f3)[0]

    if len(bad)>0:
        raise ValueError('something wrong in sampling Kroupa')
            
    x[first]=(x_min**(1.-alpha1)+ran[first]*(1.-alpha1)/k1)**(1./(1.-alpha1))
    x[second]=(x0_2**(1.-alpha2)+(1.-alpha2)/k2*(ran[second]-f1))**(1./(1.-alpha2))
    x[third]=(x0_3**(1.-alpha3)+(1.-alpha3)/k3*(ran[third]-f1-f2))**(1./(1.-alpha3))

    return x,k1,k2,k3

def bpl(size,x_min,x_max,alpha1,alpha2,x0):
    ran=np.random.uniform(low=0.,high=1.,size=size)
    x=np.zeros(size,dtype=float)

    #get normalization constant for each of three pieces
    k2_over_k1=x0**(alpha2-alpha1)
        
    if x_min<x0:
            
        if x0<=x_max:
                
            piece1=(x0**(1.-alpha1)-x_min**(1.-alpha1))/(1.-alpha1)
            piece2=(x_max**(1.-alpha2)-x0**(1.-alpha2))/(1.-alpha2)

            x0_2=x0
                                
        if x0>x_max:
                
            piece1=(x_max**(1.-alpha1)-x_min**(1.-alpha1))/(1.-alpha1)
            piece2=0.
                
            x0_2=x0
                
    if x0<=x_min:
            
        if x0<=x_max:
            piece1=0.
            piece2=(x_max**(1.-alpha2)-x_min**(1.-alpha2))/(1.-alpha2)
                
            x0_2=x_min
                
    k1=1./(piece1+piece2*k2_over_k1)#sample size normalized to 1
    k2=k1*k2_over_k1

    #get fraction of sample within each piece
    f1=k1*piece1
    f2=k2*piece2

    first=np.where(ran<f1)[0]
    second=np.where((ran>=f1)&(ran<f1+f2))[0]
    bad=np.where(ran>f1+f2)[0]

    if len(bad)>0:
        raise ValueError('something wrong in sampling BPL')
                        
    x[first]=(x_min**(1.-alpha1)+ran[first]*(1.-alpha1)/k1)**(1./(1.-alpha1))
    x[second]=(x0_2**(1.-alpha2)+(1.-alpha2)/k2*(ran[second]-f1))**(1./(1.-alpha2))
    
    return x,k1,k2

def normal_truncated(size,x_min,x_max,loc,scale):
    ran=np.random.uniform(size=size,low=0.,high=1.)
    return loc-np.sqrt(2.*scale**2)*scipy.special.erfinv(scipy.special.erf((loc-x_min)/np.sqrt(2.*scale**2))-ran*(scipy.special.erf((loc-x_min)/np.sqrt(2.*scale**2))-scipy.special.erf((loc-x_max)/np.sqrt(2.*scale**2))))

