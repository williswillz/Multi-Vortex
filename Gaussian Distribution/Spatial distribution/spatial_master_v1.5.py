import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack
import time
import matplotlib.mlab as mlab
from matplotlib.mlab import psd
from matplotlib.lines import Line2D


start = time.time()

def induced_velocity_single(x, xvort, gam, core_radius):
        """Compute velocity induced at points x by a single vortex

        Parameters
        ----------
        x : 2d array
            Locations at which to compute induced velocity.  Expressed as
            column vectors (i.e., shape should be (n,2))
        xvort : 1d array
            Location of vortex (shape should be (2,))
        gam : float
            Strength of vortex

        Notes
        -----
        Induced velocity is

        .. math:: u_\theta = \frac{\Gamma}{2 \pi r}

        where r is the distance between the point and the vortex.  If this
        distance is less than :class:`core_radius` :math:`r_0`, the velocity is
        regularized as solid-body rotation, with

        .. math:: u_\theta = \frac{\Gamma r}{2\pi r_0^2}
        """
        r = np.array(x, ndmin=2) - np.array(xvort)
        rsq = np.maximum(np.sum(r * r, 1), core_radius**2)
        # alternative regularization (Krasny, Eldredge)
        # rsq = np.sum(r * r, 1) + core_radius**2
        vel = np.transpose(np.array([-r[:,1], r[:,0]]))
        vel = gam / (2 * np.pi) * vel / rsq[:,np.newaxis]
        return np.squeeze(vel)
    
total_num = 500

xvort1 = np.zeros((total_num,2))
vel1 = np.zeros((total_num,2))
a = 0

q = np.zeros((total_num,2))
q[:,0] = np.random.uniform(0,50,total_num)
q[:,1] = np.random.uniform(-0.5,0.5,total_num)
# q[:,1] = np.random.normal(0,1,total_num)

core_radius = np.zeros((total_num,1))
core_radius[:,0] = np.random.uniform(0,1,total_num)

gam = np.zeros((total_num,1))
gam[:,0] = np.random.uniform(-1,1,total_num)

sum = np.zeros((total_num,2))

for i in range (0, total_num):
    i = i * 0.1
    xvort1[a,0] = i
    xvort1[a,1] = 0
    vel1[:,:] = induced_velocity_single(q, xvort1, gam, core_radius[:,0])
    sum[a,:] = vel1.sum(axis=0)
#    print sum 
    # vel1[:,:] = induced_velocity_single(q, xvort1, gam, core_radius)
#    print vel1
    a = a + 1

t = np.linspace(0,50,total_num)
vel_tot_mag = np.zeros((total_num,1))
vel_tot_mag = (sum[:,0]**2 + sum[:,1]**2)**0.5
#print vel_tot_mag

end = time.time()
print 'It took just '+ str(end - start) + ' seconds!'

def L_p(x):
	return 20*np.log10(np.abs(x)/2.e-5)

blocksize=100
j=0

(werte1,freq1)=psd(vel_tot_mag[:], NFFT=blocksize, Fs=100, detrend=mlab.detrend_none,window=mlab.window_hanning, noverlap=4, pad_to=None,sides='default', scale_by_freq='True')
(werte2,freq2)=psd(sum[:,0], NFFT=blocksize, Fs=100, detrend=mlab.detrend_none,window=mlab.window_hanning, noverlap=4, pad_to=None,sides='default', scale_by_freq='True')
(werte3,freq3)=psd(sum[:,1], NFFT=blocksize, Fs=100, detrend=mlab.detrend_none,window=mlab.window_hanning, noverlap=4, pad_to=None,sides='default', scale_by_freq='True')

pegel1=L_p(werte1)
pegel2=L_p(werte2)
pegel3=L_p(werte3)


fig = plt.figure(figsize=(8,40))
col = np.zeros((total_num))
col = np.array(np.where(gam<0,'r','g'))
#print col
ax0 = fig.add_subplot(611).set_xlabel('X')
ax0 = fig.add_subplot(611).set_ylabel('Y')
ax0 = fig.add_subplot(611).set_title('Vortices injection from X=0 using Normal distribution (not to scale, vortices size is 25 times the original)')
ax0 = plt.scatter(q[:,0],q[:,1], c=col[:,0],s=25*core_radius, label='vortex')
legend_elements = [Line2D([0], [0], marker='o', color='w', label='+ Vortex',
                          markerfacecolor='g', markersize=7),
                   Line2D([0], [0], marker='o', color='w', label='- Vortex',
                          markerfacecolor='r', markersize=7)]
ax0 = fig.add_subplot(611).legend(handles=legend_elements,loc='upper right')

ax1 = fig.add_subplot(612)
ax1 = fig.add_subplot(612).set_xlabel('X')
ax1 = fig.add_subplot(612).set_ylabel('Y')
ax1 = fig.add_subplot(612).set_title('Histogram for vortices injection from X=0')
ax1 = plt.hist(q[:,1], orientation='horizontal')

ax2 = fig.add_subplot(613)
ax2 = fig.add_subplot(613).set_xlabel('time')
ax2 = fig.add_subplot(613).set_ylabel('induced velocity')
ax2 = fig.add_subplot(613).set_title('induced velocity with respect to time')
ax2 = plt.plot(t, sum[:,0], c='g', label='u (x-component of induced velocity)')
ax2 = plt.plot(t, sum[:,1], c='r', label='v (y-component of induced velocity)')
#ax2 = plt.plot(t, vel_tot_mag[:], c='k', label='U (total induced velocity)')
ax2 = fig.add_subplot(613).legend(loc='best')

ax4 = fig.add_subplot(615)
ax4 = fig.add_subplot(615).set_ylabel('Power Spectral Density (dB)')
ax4 = fig.add_subplot(615).set_xlabel('Frequency (Hz)')
ax4 = fig.add_subplot(615).set_title('Power Spectral Density (PSD) of induced velocity')
ax4 = plt.semilogx(freq1,pegel1,linewidth=0.6)

ax5 = fig.add_subplot(616)
ax5 = fig.add_subplot(616).set_xlabel('X')
ax5 = fig.add_subplot(616).set_ylabel('Y')
ax5 = fig.add_subplot(616).set_title('Histogram for radius size of the vortices injected from X=0')
ax5 = plt.hist(core_radius[:], orientation='vertical')
#fig.savefig('spatial_master_v1.2.pdf')
plt.show()

fig3 = plt.figure(figsize=(15,40))
col = np.zeros((total_num))
col = np.array(np.where(gam<0,'r','g'))
#print col
ax0 = fig3.add_subplot(311).set_xlabel('X')
ax0 = fig3.add_subplot(311).set_ylabel('Y')
ax0 = fig3.add_subplot(311).set_title('Vortices injection from X=0 using Normal distribution (not to scale, vortices size is 25 times the original)')
ax0 = plt.scatter(q[:,0],q[:,1], c=col[:,0],s=25*core_radius, label='vortex')
legend_elements = [Line2D([0], [0], marker='o', color='w', label='+ Vortex',
                          markerfacecolor='g', markersize=7),
                   Line2D([0], [0], marker='o', color='w', label='- Vortex',
                          markerfacecolor='r', markersize=7)]
ax0 = fig3.add_subplot(311).legend(handles=legend_elements,loc='upper right')

ax2 = fig3.add_subplot(312)
ax2 = fig3.add_subplot(312).set_xlabel('time')
ax2 = fig3.add_subplot(312).set_ylabel('induced velocity')
ax2 = fig3.add_subplot(312).set_title('induced velocity with respect to time')
ax2 = plt.plot(t, sum[:,0], c='g', label='u (x-component of induced velocity)')
ax2 = plt.plot(t, sum[:,1], c='r', label='v (y-component of induced velocity)')
#ax2 = plt.plot(t, vel_tot_mag[:], c='k', label='U (total induced velocity)')
ax2 = fig3.add_subplot(312).legend(loc='best')


ax1 = fig3.add_subplot(313)
ax1 = fig3.add_subplot(313).set_ylabel('Power Spectral Density (dB)')
ax1 = fig3.add_subplot(313).set_xlabel('Frequency (Hz)')
ax1 = fig3.add_subplot(313).set_title('Power Spectral Density (PSD) of induced velocity')
ax1 = plt.semilogx(freq1,pegel1,c='k', label='U_total')
ax1 = plt.semilogx(freq2,pegel2,c='g', label='u_x')
ax1 = plt.semilogx(freq3,pegel3,c='r', label='v_y')
ax1 = fig3.add_subplot(313).legend(loc='best')
fig3.savefig('corrected.pdf')
plt.show()