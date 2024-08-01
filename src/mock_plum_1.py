import numpy as np
import tkinter
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy
import mycode as mc
from astropy.io import fits
import scipy.optimize
from scipy.spatial import distance
from scipy.integrate import quad
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy import interpolate
from astropy.wcs import WCS
from os import path
import scipy.ndimage.filters
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from shapely.geometry import Point,Polygon,MultiPoint
import xymass
import dill as pickle

###toy for learning 3d grid
grid_y,grid_x=np.meshgrid(np.linspace(-1.,1.,5),np.linspace(-1.,1.,4))
grid3d_y,grid3d_x,grid3d_logs=np.meshgrid(np.linspace(-1.,1.,5),np.linspace(-1.,1.,4),np.linspace(np.log10(0.),np.log10(1.),3))
grid_r=np.sqrt(grid_x**2+grid_y**2)
grid3d_r=np.sqrt(grid3d_x**2+grid3d_y**2)

grid_mask=np.zeros((4,5),dtype='bool')
grid_mask[grid_r>1.]=True

#how to apply the same mask to grid3d

#np.pause()

get_sample=True

directory='/hildafs/projects/phy200028p/mgwalker/dra2_jwst/'
rmax=3.#arcmin, matched to be slightly larger than JWST/NIRCAM FOV
mask_file=directory+'dra2_jwst_mask.pkl'
log_integral=False#perform integral in log space? (set to False for linear space)
sb_model='plum'
sep2d_model='bpl'

sep_bins=45

nplum=1000
nuni=0.*nplum
model='plum'
r_scale_true=3.*u.arcmin
dx_true=0.
dy_true=-0.
ellipticity_true=0.
position_angle_true=30.
bigsigma0_mem_true=nplum/np.pi/r_scale_true**2/(1.-ellipticity_true)

alpha_true=np.nan
beta_true=np.nan
gamma_true=np.nan

s_binary_min=0.05/60.#(0.05*u.arcsec).to(u.arcmin)#minimum separation to define binary fraction (arcmin) 
s_binary_max=10./60.#(10.*u.arcsec).to(u.arcmin)#maximum separation to define binary fraction (arcmin)
s_obs_min=0.05/60.#(0.05*u.arcsec).to(u.arcmin)#minimum separation to define observed separation function (arcmin)
s_obs_max=10./60.#(10.*u.arcsec).to(u.arcmin)#maximum separation to define observed separation function (arcmin)
s_break_1_true=2.05/60.#(2.05*u.arcsec).to(u.arcmin)
s_break_2_true=8./60.#(8.*u.arcsec).to(u.arcmin)
alpha1_1_true=2.01
alpha2_1_true=5.5
alpha1_2_true=1.1
alpha2_2_true=4.2

fb_1_true=0.5#08
fb_2_true=0.#01

bigsigma0_mem_null=nplum*(1.+fb_1_true)/np.pi/r_scale_true**2/(1.-ellipticity_true)

bigsigma0_non_true=nuni/(np.pi*(rmax*np.sqrt(3.))**2)

if get_sample:

    m_min=0.1#Msun
    m_max=np.inf#Msun
    
    plum2d=xymass.sample_r2d(size=int(nplum),model=model,ellipticity=ellipticity_true,position_angle=position_angle_true,r_scale=r_scale_true)
    uni2d=xymass.sample_r2d(size=int(nuni),model='uni',r_scale=rmax*np.sqrt(3.))

    plum_mass_primary=xymass.sample_imf(size=nplum,model='kroupa',m_min=m_min,m_max=m_max).mass*u.M_sun
    plum_mass_secondary=xymass.sample_imf(size=nplum,model='kroupa',m_min=m_min,m_max=m_max).mass*u.M_sun
    
    plum2d_with_binaries=xymass.add_binaries_func(plum2d.r_xyz,separation_func='bpl',mass_primary=plum_mass_primary,mass_secondary=plum_mass_secondary,f_binary=fb_1_true,s_min=s_obs_min,s_max=s_obs_max,alpha1=alpha1_1_true,alpha2=alpha2_1_true,s_break=s_break_1_true,projected=True)

    np.pause()
    
    
    plum_s_bpl=mc.sample_separation(int(nplum*fb_1_true),s_min=s_binary_min,s_max=s_binary_max,model='bpl',s1=s_break_1_true,alpha1=alpha1_1_true,alpha2=alpha2_1_true)#units of separations are arcmin
    plum_s_opik=mc.sample_separation(int(nplum*fb_1_true),s_min=s_binary_min,s_max=s_binary_max,model='opik')#units of separations are arcmin
    uni_s_bpl=mc.sample_separation(int(nuni*fb_2_true),s_min=s_binary_min,s_max=s_binary_max,model='bpl',s1=s_break_2_true,alpha1=alpha1_2_true,alpha2=alpha2_2_true)#units of separations are arcmin

    plum_binary_x=plum2d.x[0:int(nplum*fb_1_true)]+dx_true+plum_s_bpl.x2d
    plum_binary_y=plum2d.y[0:int(nplum*fb_1_true)]+dy_true+plum_s_bpl.y2d
    uni_binary_x=uni2d.x[0:int(nuni*fb_2_true)]+uni_s_bpl.x2d
    uni_binary_y=uni2d.y[0:int(nuni*fb_2_true)]+uni_s_bpl.y2d

    x0=np.concatenate((plum2d.x+dx_true,uni2d.x,plum_binary_x,uni_binary_x))
    y0=np.concatenate((plum2d.y+dy_true,uni2d.y,plum_binary_y,uni_binary_y))
    pop0=np.concatenate((np.full(len(plum2d.x),1),np.full(len(uni2d.x),2),np.full(len(plum_binary_x),3),np.full(len(uni_binary_x),4)))

    plum_keep=np.where((np.abs(plum2d.x)<=rmax)&(np.abs(plum2d.y)<=rmax))[0]
    uni_keep=np.where((np.abs(uni2d.x)<=rmax)&(np.abs(uni2d.y)<=rmax))[0]
    plum_binary_keep=np.where((np.abs(plum_binary_x)<=rmax)&(np.abs(plum_binary_y)<=rmax))[0]
    uni_binary_keep=np.where((np.abs(uni_binary_x)<=rmax)&(np.abs(uni_binary_y)<=rmax))[0]
    keep0=np.where((np.abs(x0)<=rmax)&(np.abs(y0)<=rmax))[0]

    x=x0[keep0]
    y=y0[keep0]
    pop=pop0[keep0]

    #perm=np.random.permutation(len(pop))
    #x=x[perm]
    #y=y[perm]
    #pop=pop[perm]

    pickle.dump((x,y,pop),open(directory+'mock_plum_1.pkl','wb'))
x,y,pop=pickle.load(open(directory+'mock_plum_1.pkl','rb'))

gridfrac_file=directory+'dra2_jwst_gridfrac.pkl'
gridfrac_logsep2d,grid_frac0,grid_frac,grid_frac2,grid_frac3=pickle.load(open(gridfrac_file,'rb'))
#grid_frac0=np.ones(np.shape(grid_frac[0]))#ones array corresponding to grid_frac when s=0
#grid_frac=np.ones(np.shape(grid_frac))#ones array corresponding to grid_frac when s=0
#grid_frac2=np.ones(np.shape(grid_frac))#ones array corresponding to grid_frac when s=0
#grid_frac3=np.ones(np.shape(grid_frac))#ones array corresponding to grid_frac when s=0
gridfrac_logsep2d=np.linspace(s_obs_min,s_obs_max,len(grid_frac))#grid of separations to use for pre-calculating completeness fractions g(R)

npix=101#np.max([1001,int(2.*rmax/0.1)+1])
nx,ny=(npix,npix)
grid_y,grid_x=np.meshgrid(np.linspace(-rmax,rmax,nx),np.linspace(-rmax,rmax,ny))
grid3d_y,grid3d_x,grid3d_logs=np.meshgrid(np.linspace(-rmax,rmax,nx),np.linspace(-rmax,rmax,ny),np.linspace(np.log10(s_obs_min),np.log10(s_obs_max),sep_bins))

grid_mask=np.full((npix,npix),False,dtype='bool')
grid3d_mask=np.concatenate(([np.expand_dims(grid_mask.T,0) for q in range(0,np.shape(grid3d_logs)[2])])).T
#crap,grid_mask=pickle.load(open(mask_file,'rb'))#use grid mask from real dra2 observation
dgrid2=(np.abs(grid_y[0][1]-grid_y[0][0]))*(np.abs(grid_x[0][0]-grid_x[1][0]))
dgrid3=(np.abs(grid_y[0][1]-grid_y[0][0]))*(np.abs(grid_x[0][0]-grid_x[1][0])*(np.abs(10.**grid3d_logs[0][0][1]-10.**grid3d_logs[0][0][0])))
dsgrid=np.abs(10.**grid3d_logs[0][0][1]-10.**grid3d_logs[0][0][0])

xxx=scipy.stats.binned_statistic_2d(x,y,x,'count',bins=npix,range=[[-rmax,rmax],[-rmax,rmax]],expand_binnumbers=True)
obs_binnumber=xxx.binnumber.T
obs_mask=mc.get_obs_mask_2dimage(obs_binnumber,grid_mask)
keep=np.where(obs_mask==False)[0]

fmemfield_true=len(np.where(pop[keep]==1)[0])/(len(np.where(pop[keep]==1)[0])+len(np.where(pop[keep]==2)[0]))
sb_prefix='/hildafs/projects/phy200028p/mgwalker/dra2_jwst/chains/sb_mock_plum_1'

sb_prior=np.array([[-2.,2.],[-2.,2.],[-10.,10.],[0.,0.99999999],[-90.,90.],[-1.,3.],[0.,0.]])#(dx,dy,log10_bigsigma0_mem,fmem,PA,log10_r_scale,ellipticity)
if ((model=='plum_n')|(model=='exp_n')):
    sb_prior=np.array([[-2.,2.],[-2.,2.],[-10.,10.],[0.,0.99999999],[-90.,90.],[-1.,3.],[0.,0.],[0.,1.]])#(dx,dy,log10_bigsigma0_mem,fmem,PA,log10_r_scale,ellipticity,n)
if (model=='2bg'):
    sb_prior=np.array([[-2.,2.],[-2.,2.],[-10.,10.],[0.,0.99999999],[-90.,90.],[-1.,3.],[0.,0.],[3.1,15.],[0.,1.999]])#(dx,dy,log10_bigsigma0_mem,fmem,PA,log10_r_scale,ellipticity,beta,gamma)
if (model=='abg'):
    sb_prior=np.array([[-2.,2.],[-2.,2.],[-10.,10.],[0.,0.99999999],[-90.,90.],[-1.,3.],[0.,0.],[0.51,3.],[3.1,15.],[0.,2.]])#(dx,dy,log10_bigsigma0_mem,fmem,PA,log10_r_scale,ellipticity,alpha,beta,gamma)
    
sb_resume=True
sb_result,sb_bestfit=mc.fit_sb_grid(x[keep],y[keep],grid_x,grid_y,prior=sb_prior,obs_mask=obs_mask[keep],grid_mask=grid_mask,model=model,prefix=sb_prefix,resume=sb_resume,interp=False,sampler='multinest')

#sb_grid_pdfname='poo.pdf'
#mc.get_sb_grid_figure(x[keep],y[keep],grid_x,grid_y,obs_mask[keep],grid_mask,sb_model,sb_bestfit,rmax,nx,r_scale_true,ellipticity_true,position_angle_true,sb_grid_pdfname,show=True)

sb_grid_true,sb_grid_mem_true,crap,nmemfield_true,bigsigma0_non_true=mc.get_sb_model_grid(grid_x,grid_y,obs_mask=obs_mask,grid_mask=grid_mask,model=model,dx=dx_true,dy=dy_true,bigsigma0_mem=bigsigma0_mem_true,fmemfield=fmemfield_true,ellipticity=ellipticity_true,position_angle=position_angle_true,r_scale=r_scale_true,n_index=np.nan,alpha=np.nan,beta=np.nan,gamma=np.nan)
sb_grid_null,sb_grid_mem_null,crap,nmemfield_null,bigsigma0_non_null=mc.get_sb_model_grid(grid_x,grid_y,obs_mask=obs_mask,grid_mask=grid_mask,model=model,dx=dx_true,dy=dy_true,bigsigma0_mem=bigsigma0_mem_null,fmemfield=fmemfield_true,ellipticity=ellipticity_true,position_angle=position_angle_true,r_scale=r_scale_true,n_index=np.nan,alpha=np.nan,beta=np.nan,gamma=np.nan)

bigsigma0_non_true=nuni/(np.pi*(rmax*np.sqrt(3.))**2)

sb_grid_non_true=sb_grid_true-sb_grid_mem_true
nnonfield_true=dgrid2*np.sum(sb_grid_non_true[grid_mask==False])

sb_grid_non_null=sb_grid_null-sb_grid_mem_null
nnonfield_null=dgrid2*np.sum(sb_grid_non_null[grid_mask==False])

sample_xy0=np.c_[x[keep],y[keep]]
sample_tree=scipy.spatial.KDTree(sample_xy0)

logs_binary=np.linspace(np.log10(s_binary_min),np.log10(s_binary_max),sep_bins+1)
logs_obs=np.linspace(np.log10(s_obs_min),np.log10(s_obs_max),sep_bins+1)
s_binary=10.**logs_binary
s_obs=10.**logs_obs

phib_1_true=mc.sample_separation(size=1,model='bpl',alpha1=alpha1_1_true,alpha2=alpha2_1_true,s1=s_break_1_true,s_min=s_binary_min,s_max=s_binary_max)
phib_2_true=mc.sample_separation(size=1,model='bpl',alpha1=alpha1_2_true,alpha2=alpha2_2_true,s1=s_break_2_true,s_min=s_binary_min,s_max=s_binary_max)

nn,nn_partner,sep2d,nn_hist,sep2d_hist,lognn_hist,logsep2d_hist,pair_obs_xy=mc.get_separations(x[keep],y[keep],s_min=s_obs_min,s_max=s_obs_max,bins=sep_bins)

pair_x=[pair_obs_xy[q][0][0] for q in range(0,len(pair_obs_xy))]
pair_y=[pair_obs_xy[q][0][1] for q in range(0,len(pair_obs_xy))]
sample_xys=np.c_[pair_x,pair_y,np.log10(sep2d)]

bins1=np.linspace(-rmax,rmax,npix+1)
bins2=np.linspace(np.log10(s_obs_min),np.log10(s_obs_max),sep_bins+1)

pair_obs_keep=np.where((sep2d>=s_obs_min)&(sep2d<=s_obs_max))[0]
xxx=scipy.stats.binned_statistic_dd(sample_xys[pair_obs_keep],sep2d[pair_obs_keep],statistic='count',bins=[bins1,bins1,bins2],expand_binnumbers=True)
yyy=scipy.stats.binned_statistic_2d(x[keep],y[keep],x[keep],'count',bins=[bins1,bins1],expand_binnumbers=True)

phib1=phib_1_true.func(10.**grid3d_logs.flatten()).reshape(npix,npix,sep_bins)
phib2=phib_2_true.func(10.**grid3d_logs.flatten()).reshape(npix,npix,sep_bins)

gridfrac_interp=1.

grid3d_bigsigma,grid3d_bigsigma_mem,e_grid3d_bigsigma,nmemfield,bigsigma0_non=mc.get_sb_model_grid3d(grid3d_x,grid3d_y,grid3d_mask=grid3d_mask,model=sb_model,dx=dx_true,dy=dy_true,bigsigma0_mem=bigsigma0_mem_true,bigsigma0_non=bigsigma0_non_true,ellipticity=ellipticity_true,position_angle=-position_angle_true,r_scale=r_scale_true,n_index=0.)

fbtry=np.linspace(0.,1.,100)
loglike=[]
for i in range(0,len(fbtry)):
    grid3d_bigsigma,grid3d_bigsigma_mem,e_grid3d_bigsigma,nmemfield,bigsigma0_non=mc.get_sb_model_grid3d(grid3d_x,grid3d_y,grid3d_mask=grid3d_mask,model=sb_model,dx=dx_true,dy=dy_true,bigsigma0_mem=bigsigma0_mem_true,bigsigma0_non=bigsigma0_non_true,ellipticity=ellipticity_true,position_angle=position_angle_true,r_scale=r_scale_true,n_index=0.)
    term1=2.*fbtry[i]*grid3d_bigsigma_mem*phib1+2.*fb_2_true*bigsigma0_non*phib2
    term2=(1.+fbtry[i]+fbtry[i]+fbtry[i]*fbtry[i])*2.*np.pi*10.**grid3d_logs*grid3d_bigsigma_mem*grid3d_bigsigma_mem*gridfrac_interp
    term3=(1.+fbtry[i]+fb_2_true+fbtry[i]*fb_2_true)*2.*np.pi*10.**grid3d_logs*grid3d_bigsigma_mem*bigsigma0_non*gridfrac_interp
    term4=(1.+fb_2_true+fbtry[i]+fb_2_true*fbtry[i])*2.*np.pi*10.**grid3d_logs*bigsigma0_non*grid3d_bigsigma_mem*gridfrac_interp
    term5=(1.+fb_2_true+fb_2_true+fb_2_true*fb_2_true)*2.*np.pi*10.**grid3d_logs*bigsigma0_non*bigsigma0_non*gridfrac_interp
    grid3d_n=(term1+term2+term3+term4+term5)*dgrid3
    nnn=np.sqrt(grid3d_n)
    nnnn=(nnn*nnn)/2.
    pois=scipy.stats.poisson.pmf(xxx.statistic.flatten(),nnnn.flatten())
    loglike.append(np.sum(scipy.stats.poisson.logpmf(xxx.statistic[grid3d_mask==False].flatten(),nnnn[grid3d_mask==False].flatten())))
    print(fbtry[i],loglike[len(loglike)-1])
loglike=np.array(loglike)
plt.plot(fbtry,loglike)
plt.show()

grid3d_bigsigma,grid3d_bigsigma_mem,e_grid3d_bigsigma,nmemfield,bigsigma0_non=mc.get_sb_model_grid3d(grid3d_x,grid3d_y,grid3d_mask=grid3d_mask,model=sb_model,dx=dx_true,dy=dy_true,bigsigma0_mem=bigsigma0_mem_true,bigsigma0_non=bigsigma0_non_true,ellipticity=ellipticity_true,position_angle=-position_angle_true,r_scale=r_scale_true,n_index=0.)
term1=2.*fb_1_true*grid3d_bigsigma_mem*phib1+2.*fb_2_true*bigsigma0_non*phib2
term2=(1.+fb_1_true+fb_1_true+fb_1_true*fb_1_true)*2.*np.pi*10.**grid3d_logs*grid3d_bigsigma_mem*grid3d_bigsigma_mem*gridfrac_interp
term3=(1.+fb_1_true+fb_2_true+fb_1_true*fb_2_true)*2.*np.pi*10.**grid3d_logs*grid3d_bigsigma_mem*bigsigma0_non*gridfrac_interp
term4=(1.+fb_2_true+fb_1_true+fb_2_true*fb_1_true)*2.*np.pi*10.**grid3d_logs*bigsigma0_non*grid3d_bigsigma_mem*gridfrac_interp
term5=(1.+fb_2_true+fb_2_true+fb_2_true*fb_2_true)*2.*np.pi*10.**grid3d_logs*bigsigma0_non*bigsigma0_non*gridfrac_interp
grid3d_n=(term1+term2+term3+term4+term5)*dgrid3
pois=scipy.stats.poisson.pmf(xxx.statistic.flatten(),grid3d_n.flatten())
loglike=np.sum(scipy.stats.poisson.logpmf(xxx.statistic[grid3d_mask==False].flatten(),grid3d_n[grid3d_mask==False].flatten()))

grid_npair=np.sum(grid3d_n,axis=2)
grid_n=sb_grid_true*dgrid2

cube_true=np.array([fb_1_true,fb_2_true,np.log10(s_break_1_true),np.log10(s_break_2_true),alpha1_1_true,alpha1_2_true,alpha2_1_true-alpha1_1_true,alpha2_2_true-alpha1_2_true,dx_true,dy_true,np.log10(bigsigma0_mem_true*(1.+fb_1_true)),np.log10(bigsigma0_non_true*(1.+fb_2_true)),position_angle_true,np.log10(r_scale_true),ellipticity_true])
loglike_true=mc.fit_separations_loglike_fudge0(cube_true,grid_x,grid_y,obs_mask,grid_mask,model,s_binary_min,s_binary_max,sep2d_model,sb_model,logsep2d_hist[1],logsep2d_hist[0],dgrid2,gridfrac_logsep2d,grid_frac,grid_frac0,grid_frac2,grid_frac3,use_fmemfield=False)

fb=np.linspace(0.,1.,sep_bins)
loglike=[]
loglike_fudge=[]
for i in range(0,len(fb)):
    print(i)
    cube=np.array([fb[i],fb_2_true,np.log10(s_break_1_true),np.log10(s_break_2_true),alpha1_1_true,alpha1_2_true,alpha2_1_true-alpha1_1_true,alpha2_2_true-alpha1_2_true,dx_true,dy_true,np.log10(bigsigma0_mem_true*(1.+fb[i])),np.log10(bigsigma0_non_true*(1.+fb_2_true)),position_angle_true,np.log10(r_scale_true),ellipticity_true])
    #cube=np.array([fb[i],fb_2_true,np.log10(s_break_1_true),np.log10(s_break_2_true),alpha1_1_true,alpha1_2_true,alpha2_1_true-alpha1_1_true,alpha2_2_true-alpha1_2_true,dx_true,dy_true,np.log10(bigsigma0_mem_true),np.log10(bigsigma0_non_true),position_angle_true,np.log10(r_scale_true),ellipticity_true])
    if log_integral:
        loglike.append(mc.fit_separations_loglike0(cube,grid_x,grid_y,obs_mask,grid_mask,model,s_binary_min,s_binary_max,sep2d_model,sb_model,logsep2d_hist[1],logsep2d_hist[0],dgrid2,gridfrac_logsep2d,grid_frac,grid_frac0,grid_frac2,grid_frac3,log_integral=True))
    else:
        #loglike.append(mc.fit_separations_loglike0(cube,grid_x,grid_y,obs_mask,grid_mask,model,s_binary_min,s_binary_max,sep2d_model,sb_model,logsep2d_hist[1],logsep2d_hist[0],dgrid2,gridfrac_logsep2d,grid_frac,grid_frac0,grid_frac2,grid_frac3,use_fmemfield=False))
        loglike_fudge.append(mc.fit_separations_loglike_fudge0(cube,grid_x,grid_y,obs_mask,grid_mask,model,s_binary_min,s_binary_max,sep2d_model,sb_model,logsep2d_hist[1],logsep2d_hist[0],dgrid2,gridfrac_logsep2d,grid_frac,grid_frac0,grid_frac2,grid_frac3,use_fmemfield=False))

loglike=np.array(loglike)
loglike_fudge=np.array(loglike_fudge)

#plt.plot(fb,loglike,label='more exact')
plt.plot(fb,loglike_fudge,label='fudge')
plt.axvline(fb_1_true,linestyle=':',color='k')
plt.legend(loc=4)
plt.show()
plt.close()

sb_bigsigma0_non=[]
for j in range(0,len(sb_result['samples'])):
    dx=sb_result['samples'][j][0]
    dy=sb_result['samples'][j][1]
    #bigsigma0_mem=(1.+sb_result['samples'][j][0])*10.**sb_result['samples'][j][2]
    bigsigma0_mem=10.**sb_result['samples'][j][2]
    fmemfield=sb_result['samples'][j][3]
    ellipticity=sb_result['samples'][j][6]
    position_angle=sb_result['samples'][j][4]
    r_scale=10.**sb_result['samples'][j][5]        

    sb_grid,sb_grid_mem,crap,nmemfield,bigsigma0_non=mc.get_sb_model_grid(grid_x,grid_y,obs_mask=obs_mask,grid_mask=grid_mask,model=model,dx=dx,dy=dy,bigsigma0_mem=bigsigma0_mem_true,fmemfield=fmemfield,ellipticity=ellipticity,position_angle=position_angle,r_scale=r_scale,n_index=0.,alpha=np.nan,beta=np.nan,gamma=np.nan)
    sb_bigsigma0_non.append(bigsigma0_non)
sb_bigsigma0_non=np.array(sb_bigsigma0_non)

separations_prior=np.array([[0.,1.],[0.,1.],[np.log10(s_binary_min),np.log10(s_binary_max)],[np.log10(s_binary_min),np.log10(s_binary_max)],[1.001,5.],[1.001,5.],[0.,5.],[0.,5.],[np.mean(sb_result['samples'].T[0]),np.std(sb_result['samples'].T[0])],[np.mean(sb_result['samples'].T[1]),np.std(sb_result['samples'].T[1])],[np.mean(sb_result['samples'].T[2]),np.std(sb_result['samples'].T[2])],[np.mean(np.log10(sb_bigsigma0_non)),np.std(np.log10(sb_bigsigma0_non))],[np.mean(sb_result['samples'].T[4]),np.std(sb_result['samples'].T[4])],[np.mean(sb_result['samples'].T[5]),np.std(sb_result['samples'].T[5])],[np.mean(sb_result['samples'].T[6]),np.std(sb_result['samples'].T[6])]])#fb_1,fb_2,log10[sbreak1/arcmin],log10[sbreak2/arcmin],alph1_1,alpha1_2,alpha2_1,alpha2_2,dx,dy,log10bigsigma0_mem',log10bigsigma0_non',ellipticity,position_angle,log10r_scale
separations_prefix=directory+'chains/separations_mock_plum_1'
separations_resume=True
np.pause()
separations_result,separations_bestfit=mc.fit_separations_grid(logsep2d_hist[1],logsep2d_hist[0],s_binary_min,s_binary_max,sep2d_model,sb_model,sb_result,grid_x,grid_y,grid_mask,obs_mask,gridfrac_logsep2d,grid_frac0,grid_frac,grid_frac2,grid_frac3,prior=separations_prior,prefix=separations_prefix,resume=separations_resume)

fb_1_bestfit=separations_bestfit['parameters'][0]
fb_2_bestfit=separations_bestfit['parameters'][1]
s_break_1_bestfit=10.**separations_bestfit['parameters'][2]
s_break_2_bestfit=10.**separations_bestfit['parameters'][3]
alpha1_1_bestfit=separations_bestfit['parameters'][4]
alpha1_2_bestfit=separations_bestfit['parameters'][5]
alpha2_1_bestfit=separations_bestfit['parameters'][6]+alpha1_1_bestfit
alpha2_2_bestfit=separations_bestfit['parameters'][7]+alpha1_2_bestfit
dx_bestfit=separations_bestfit['parameters'][8]
dy_bestfit=separations_bestfit['parameters'][9]
bigsigma0_mem_bestfit=10.**separations_bestfit['parameters'][10]/(1.+fb_1_bestfit)
bigsigma0_non_bestfit=10.**separations_bestfit['parameters'][11]/(1.+fb_2_bestfit)
position_angle_bestfit=separations_bestfit['parameters'][12]
r_scale_bestfit=10.**separations_bestfit['parameters'][13]
ellipticity_bestfit=separations_bestfit['parameters'][14]

sb_grid_bestfit,sb_grid_mem_bestfit,crap,nmemfield_bestfit,bigsigma0_non_bestfit=mc.get_sb_model_grid(grid_x,grid_y,obs_mask=obs_mask,grid_mask=grid_mask,model=model,dx=dx_bestfit,dy=dy_bestfit,bigsigma0_mem=bigsigma0_mem_bestfit,bigsigma0_non=bigsigma0_non_bestfit,ellipticity=ellipticity_bestfit,position_angle=position_angle_bestfit,r_scale=r_scale_bestfit,n_index=0.)
sb_grid_non_bestfit=sb_grid_bestfit-sb_grid_mem_bestfit
nnonfield_bestfit=dgrid2*np.sum(sb_grid_non_bestfit[grid_mask==False])

bestfit_cube=separations_bestfit['parameters']
bestfit2_cube=np.copy(bestfit_cube)
bestfit2_cube[8:14]=cube_true[8:14]
#bestfit_cube=np.array([fb_1_bestfit,fb_2_bestfit,np.log10(s_break_1_bestfit),np.log10(s_break_2_bestfit),alpha1_1_bestfit,alpha1_2_bestfit,alpha2_1_bestfit-alpha1_1_bestfit,alpha2_2_bestfit-alpha1_2_bestfit,dx_bestfit,dy_bestfit,np.log10(bigsigma0_mem_bestfit),np.log10(bigsigma0_non_bestfit),ellipticity_bestfit,position_angle_bestfit,np.log10(r_scale_bestfit)])

bestfit_loglike=mc.fit_separations_loglike_fudge0(bestfit_cube,grid_x,grid_y,obs_mask,grid_mask,model,s_binary_min,s_binary_max,sep2d_model,sb_model,logsep2d_hist[1],logsep2d_hist[0],dgrid2,gridfrac_logsep2d,grid_frac,grid_frac0,grid_frac2,grid_frac3,use_fmemfield=False)

phib_1_bestfit=mc.sample_separation(size=1,model='bpl',alpha1=alpha1_1_bestfit,alpha2=alpha2_1_bestfit,s1=s_break_1_bestfit,s_min=s_binary_min,s_max=s_binary_max)
phib_2_bestfit=mc.sample_separation(size=1,model='bpl',alpha1=alpha1_2_bestfit,alpha2=alpha2_2_bestfit,s1=s_break_2_bestfit,s_min=s_binary_min,s_max=s_binary_max)

ps=[]
ps0=[]
ps_bestfit=[]
ps_null=[]
for i in range(1,len(s_obs)):
    low=s_obs[i-1]
    high=s_obs[i]
    middle=(low+high)/2.
    #I=scipy.integrate.quad(mc.separations_integrand,low,high,args=(sb_grid_mem_true,sb_grid_non_true,grid_mask,dgrid2,gridfrac_logsep2d,grid_frac0,grid_frac,grid_frac2,grid_frac3,fb_1_true,fb_2_true,phib_1_true,phib_2_true))
    #n=np.sqrt(I[0])
    I=(high-low)*mc.get_psep(middle,sb_grid_mem_true,sb_grid_non_true,grid_mask,dgrid2,gridfrac_logsep2d,grid_frac0,grid_frac,grid_frac2,grid_frac3,fb_1_true,fb_2_true,phib_1_true,phib_2_true)
    n=np.sqrt(I)
    print(i,low*60,high*60,(n*n)/2.)
    ps.append((n*n)/2.)#number of distinct pairs (not involving same star) is (N^2-N)/2

    #I=scipy.integrate.quad(mc.separations_integrand,low,high,args=(sb_grid_mem_bestfit,sb_grid_non_bestfit,grid_mask,dgrid2,gridfrac_logsep2d,grid_frac0,grid_frac,grid_frac2,grid_frac3,fb_1_bestfit,fb_2_bestfit,phib_1_bestfit,phib_2_bestfit))
    #n=np.sqrt(I[0])
    I=(high-low)*mc.get_psep(middle,sb_grid_mem_bestfit,sb_grid_non_bestfit,grid_mask,dgrid2,gridfrac_logsep2d,grid_frac0,grid_frac,grid_frac2,grid_frac3,fb_1_bestfit,fb_2_bestfit,phib_1_bestfit,phib_2_bestfit)
    n=np.sqrt(I)
    ps_bestfit.append((n*n)/2.)#number of distinct pairs (not involving same star) is (N^2-N)/2
    
    #I=scipy.integrate.quad(mc.separations_integrand,low,high,args=(sb_grid_mem_null,sb_grid_non_null,grid_mask,dgrid2,gridfrac_logsep2d,grid_frac0,grid_frac,grid_frac2,grid_frac3,0.,0.,phib_1_true,phib_2_true))
    #n=np.sqrt(I[0])
    I=(high-low)*mc.get_psep(middle,sb_grid_mem_true,sb_grid_non_true,grid_mask,dgrid2,gridfrac_logsep2d,grid_frac0,grid_frac,grid_frac2,grid_frac3,0.,0.,phib_1_true,phib_2_true)
    n=np.sqrt(I)
    ps_null.append((n*n)/2.)#number of distinct pairs (not involving same star) is (N^2-N)/2
    
    #low=s_binary[i-1]
    #high=s_binary[i]
    #I=scipy.integrate.quad(mc.separations_integrand,low,high,args=(sb_grid_mem,sb_grid_non,grid_mask,dgrid2,gridfrac_logsep2d,grid_frac0,grid_frac,grid_frac2,grid_frac3,fb_1,fb_2,phib_1,phib_2))
    #n=np.sqrt(I[0])
    #ps0.append((n*n-n)/2.)#number of distinct pairs (not involving same star) is (N^2-N)/2    
ps=np.array(ps)
ps_bestfit=np.array(ps_bestfit)
#ps0=np.array(ps0)
ps_null=np.array(ps_null)

#plt.plot(s_binary[1:]*60,ps0,color='r',linestyle=':')
bincenters=np.array([(10.**logsep2d_hist[1][i-1]+10.**logsep2d_hist[1][i])/2. for i in range(1,len(logsep2d_hist[1]))])
plt.step(10.**logsep2d_hist[1][0:len(logsep2d_hist[1])-1]*60.,logsep2d_hist[0],color='k',where='post')#
#plt.step(bincenters*60.,logsep2d_hist[0],color='k',where='pre')#
#plt.show()
#plt.scatter(s_obs[1:]*60,ps,color='k',s=10)
plt.scatter(bincenters*60,ps,color='k',s=10)
plt.plot(bincenters*60,ps,color='k',linestyle=':',label='truth')
plt.plot(bincenters*60,ps_bestfit,color='r',linestyle=':',label='best fit')
plt.plot(bincenters*60,ps_null,color='b',linestyle=':',label='null')
#plt.scatter(s_binary[1:]*60,ps0,color='r',s=5)
plt.axvline(s_binary_min*60,linestyle='--',color='r')
plt.axvline(s_binary_max*60,linestyle='--',color='r')
plt.axvline(s_obs_min*60,linestyle='--',color='k')
plt.axvline(s_obs_max*60,linestyle='--',color='k')
plt.xlabel('2d separation [arcsec]')
plt.ylabel('N')
plt.yscale('log')
plt.legend(loc=2)
plt.xscale('log')
plt.yscale('log')
plt.show()
plt.close()

#pars=[]
#with open(directory+'chains/crap2dra2_jwst_separationspost_equal_weights.dat') as f:
#    data=f.readlines()
#for line in data:
#    p=line.split()
#    pars.append([float(p[0]),float(p[1]),float(p[2]),float(p[3]),float(p[4]),float(p[5]),float(p[6]),float(p[7]),float(p[8]),float(p[9]),float(p[10]),float(p[11]),float(p[12]),float(p[13]),float(p[14]),float(p[15])])
#pars=np.array(pars)

shite=separations_result['samples']
#shite=pars

ps=[]
ps_null=[]
dgrid2=(np.abs(grid_y[0][1]-grid_y[0][0]))*(np.abs(grid_x[0][0]-grid_x[1][0]))

#for j in range(0,len(sb_result['samples'])):
for j in range(1050,1060):
    dx=sb_result['samples'][j][0]
    dy=sb_result['samples'][j][1]
    #bigsigma0_mem=(1.+sb_result['samples'][j][0])*10.**sb_result['samples'][j][2]
    bigsigma0_mem=10.**sb_result['samples'][j][2]
    fmemfield=sb_result['samples'][j][3]
    ellipticity=sb_result['samples'][j][6]
    position_angle=sb_result['samples'][j][4]
    r_scale=10.**sb_result['samples'][j][5]        
    n_index=sb_result['samples'][j][7]
    sb_grid,sb_grid_mem,crap,nmemfield,bigsigma0_non=mc.get_sb_model_grid(grid_x,grid_y,obs_mask=obs_mask,grid_mask=grid_mask,model=sb_model,dx=dx,dy=dy,bigsigma0_mem=bigsigma0_mem,fmemfield=fmemfield,position_angle=position_angle,r_scale=r_scale,ellipticity=ellipticity,n_index=n_index)
    sb_grid_non=sb_grid-sb_grid_mem
    
    bigsigma0_mem_null=bigsigma0_mem#*(1.+sb_result['samples'][j][0])
    sb_grid_null,sb_grid_mem_null,crap,nmemfield_null,bigsigma0_non_null=mc.get_sb_model_grid(grid_x,grid_y,obs_mask=obs_mask,grid_mask=grid_mask,model=sb_model,dx=dx,dy=dy,bigsigma0_mem=bigsigma0_mem_null,fmemfield=fmemfield,ellipticity=ellipticity,position_angle=position_angle,r_scale=r_scale,n_index=0.)
    sb_grid_non_null=sb_grid_null-sb_grid_mem_null
    
    phib_1=mc.sample_separation(size=1,model='bpl',alpha1=np.mean(separations_prior[4]),alpha2=np.mean(separations_prior[6])+np.mean(separations_prior[4]),s1=10.**np.mean(separations_prior[2]),s_min=s_binary_min,s_max=s_binary_max)
    phib_2=mc.sample_separation(size=1,model='bpl',alpha1=np.mean(separations_prior[5]),alpha2=np.mean(separations_prior[7])+np.mean(separations_prior[5]),s1=10.**np.mean(separations_prior[3]),s_min=s_binary_min,s_max=s_binary_max)
    fb_1=0.
    fb_2=0.
        
    ps0=[]
    for i in range(1,len(logsep2d_hist[1])):
        low=10.**logsep2d_hist[1][i-1]
        high=10.**logsep2d_hist[1][i]
        middle=(low+high)/2.
        #I=scipy.integrate.quad(mc.separations_integrand,low,high,args=(sb_grid_mem_null,sb_grid_non_null,grid_mask,dgrid2,gridfrac_logsep2d,grid_frac0,grid_frac,grid_frac2,grid_frac3,fb_1,fb_2,phib_1,phib_2))
        #n=np.sqrt(I[0])
        I=(high-low)*mc.get_psep(middle,sb_grid_mem_null,sb_grid_non_null,grid_mask,dgrid2,gridfrac_logsep2d,grid_frac0,grid_frac,grid_frac2,grid_frac3,fb_1,fb_2,phib_1,phib_2)
        n=np.sqrt(I)
        ps0.append((n*n)/2.)#number of distinct pairs (not involving same star) is (N^2-N)/2

    ps0=np.array(ps0)
    ps_null.append(ps0)
ps_null=np.array(ps_null)
separations_posterior_null=ps_null

#for j in range(0,len(shite)):
for j in range(1050,1060):
    sb_grid,sb_grid_mem,crap,nmemfield,bigsigma0_non=mc.get_sb_model_grid(grid_x,grid_y,obs_mask=obs_mask,grid_mask=grid_mask,model=sb_model,dx=shite[j][8],dy=shite[j][9],bigsigma0_mem=10.**shite[j][10]/(1.+shite[j][0]),bigsigma0_non=10.**shite[j][11]/(1.+shite[j][1]),ellipticity=shite[j][14],position_angle=shite[j][12],r_scale=10.**shite[j][13],n_index=0.)
    sb_grid_non=sb_grid-sb_grid_mem
    phib_1=mc.sample_separation(size=1,model='bpl',alpha1=shite[j][4],alpha2=shite[j][6]+shite[j][4],s1=10.**shite[j][2],s_min=s_binary_min,s_max=s_binary_max)
    phib_2=mc.sample_separation(size=1,model='bpl',alpha1=shite[j][5],alpha2=shite[j][7]+shite[j][5],s1=10.**shite[j][3],s_min=s_binary_min,s_max=s_binary_max)
    fb_1=shite[j][0]
    fb_2=shite[j][1]

    ps0=[]
    for i in range(1,len(logsep2d_hist[1])):
        low=10.**logsep2d_hist[1][i-1]
        high=10.**logsep2d_hist[1][i]
        middle=(low+high)/2.
        #I=scipy.integrate.quad(mc.separations_integrand,low,high,args=(sb_grid_mem,sb_grid_non,grid_mask,dgrid2,gridfrac_logsep2d,grid_frac0,grid_frac,grid_frac2,grid_frac3,fb_1,fb_2,phib_1,phib_2))
        #n=np.sqrt(I[0])
        I=(high-low)*mc.get_psep(middle,sb_grid_mem,sb_grid_non,grid_mask,dgrid2,gridfrac_logsep2d,grid_frac0,grid_frac,grid_frac2,grid_frac3,fb_1,fb_2,phib_1,phib_2)
        n=np.sqrt(I)
        ps0.append((n*n)/2.)#number of distinct pairs (not involving same star) is (N^2-N)/2

    ps0=np.array(ps0)
    ps.append(ps0)
ps=np.array(ps)
separations_posterior=ps

shite=separations_bestfit['parameters']
#best=np.where(pars.T[15]==np.max(pars.T[15]))[0][0]
#shite=pars[best]

sb_grid,sb_grid_mem,crap,nmemfield,bigsigma0_non=mc.get_sb_model_grid(grid_x,grid_y,obs_mask=obs_mask,grid_mask=grid_mask,model=sb_model,dx=shite[8],dy=shite[9],bigsigma0_mem=10.**shite[10]/(1.+shite[0]),bigsigma0_non=10.**shite[11]/(1.+shite[1]),position_angle=shite[12],r_scale=10.**shite[13],ellipticity=shite[14],n_index=0.)
sb_grid_non=sb_grid-sb_grid_mem
phib_1=mc.sample_separation(size=1,model='bpl',alpha1=shite[4],alpha2=shite[6]+shite[4],s1=10.**shite[2],s_min=s_binary_min,s_max=s_binary_max)
phib_2=mc.sample_separation(size=1,model='bpl',alpha1=shite[5],alpha2=shite[7]+shite[5],s1=10.**shite[3],s_min=s_binary_min,s_max=s_binary_max)
fb_1=shite[0]
fb_2=shite[1]

ps_bestfit=[]
for i in range(1,len(logsep2d_hist[1])):
    low=10.**logsep2d_hist[1][i-1]
    high=10.**logsep2d_hist[1][i]
    I=scipy.integrate.quad(mc.separations_integrand,low,high,args=(sb_grid_mem,sb_grid_non,grid_mask,dgrid2,gridfrac_logsep2d,grid_frac0,grid_frac,grid_frac2,grid_frac3,fb_1,fb_2,phib_1,phib_2))
    n=np.sqrt(I[0])
    ps_bestfit.append((n*n)/2.)#number of distinct pairs (not involving same star) is (N^2-N)/2
ps_bestfit=np.array(ps_bestfit)

gs=plt.GridSpec(20,20)
gs.update(wspace=0,hspace=0)
fig=plt.figure(figsize=(6,6))
ax1=fig.add_subplot(gs[0:9,0:9])

bincenters=np.array([(10.**logsep2d_hist[1][i-1]+10.**logsep2d_hist[1][i])/2. for i in range(1,len(logsep2d_hist[1]))])

ax1.plot(bincenters*60.,[np.median(separations_posterior.T[q]) for q in range(0,len(separations_posterior.T))],color='r',alpha=0.3,label='test')
ax1.plot(bincenters*60.,ps_bestfit,color='orange',alpha=0.3,label='best fit')
ax1.fill_between(bincenters*60.,[np.percentile(separations_posterior.T[q],16) for q in range(0,len(separations_posterior.T))],[np.percentile(separations_posterior.T[q],84) for q in range(0,len(separations_posterior.T))],color='r',alpha=0.3)
ax1.plot(bincenters*60.,[np.median(separations_posterior_null.T[q]) for q in range(0,len(separations_posterior_null.T))],color='b',alpha=0.3,label='null')
ax1.fill_between(bincenters*60.,[np.percentile(separations_posterior_null.T[q],16) for q in range(0,len(separations_posterior_null.T))],[np.percentile(separations_posterior_null.T[q],84) for q in range(0,len(separations_posterior_null.T))],color='b',alpha=0.3)

ax1.step(bincenters*60,logsep2d_hist[0],color='k',where='mid')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('separation [arcsec]')
ax1.set_ylabel('number of pairs')
ax1.legend()
plt.show()
plt.close()
