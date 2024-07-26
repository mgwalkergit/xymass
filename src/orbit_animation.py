import numpy as np
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import astropy.units as u

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

