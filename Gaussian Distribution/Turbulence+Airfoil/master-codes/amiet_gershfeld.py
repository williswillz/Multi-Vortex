# -*- coding: utf-8 -*-
"""
Created on Thu May  9 19:28:10 2019

@author: WS1
"""

#Amiet formulation
#Amiet formulation with Gershfeld's correction

from pylab import *
from numpy import *

# Beschriftungseigenschaften
ts=8
rc('font', family='sans-serif',size=ts)
rc('axes', labelsize=ts)
rc('xtick',labelsize=ts)
rc('ytick',labelsize=ts)
rc('legend',fontsize=ts,handlelength=1.6,labelspacing=0.3,handletextpad=0.5,borderpad=0.3)
width=8.
height=6.
ymin,ymax=-20,120

collist=['b','g','r','c','m']



###################
##     Berechnung LEN
###################

# Amiet
def calcAmiet(bx,Rx,Ux,ux,Lambdax,tfreqx,cx):
	Mx=Ux/cx
	erg_ami=zeros(len(tfreqx))
	zaehler=0
	for fx in tfreqx:
		Kx=8.*pi*fx*Lambdax/(3.*Ux)
		erg_ami[zaehler]=10.*log10((bx*Lambdax*Mx**5.*ux**2.*Kx**3.)/(Rx**2.*Ux**2.*(1+Kx**2.)**(7./3.)))+181.3
		zaehler=zaehler+1
	return(erg_ami)

# Amiet + Gershfeld
def calcAmietGershfeld(bx,dx,Rx,Ux,ux,Lambdax,tfreqx,cx):
	Mx=Ux/cx
	erg_ami=zeros(len(tfreqx))
	zaehler=0
	for fx in tfreqx:
		Kx=8.*pi*fx*Lambdax/(3.*Ux)
		erg_ami[zaehler]=10.*log10(exp(-2.*pi*fx*dx/(2.*Ux))*(bx*Lambdax*Mx**5.*ux**2.*Kx**3.)/(Rx**2.*Ux**2.*(1+Kx**2.)**(7./3.)))+181.3
		zaehler=zaehler+1
	return (erg_ami)


#####################
##        Tu
#####################
chord=1
span=5
U=5.0
Lambda=0.0058
freq=array((10,12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,3150,4000,5000,6300,8000,10000,12500,16000,20000))

Tulist=[3.8]

labellist=[r'$Tu$ = 3.8%']


figure(2,figsize=(width/2.54,height/2.54))
for Tu in Tulist:
	j=Tulist.index(Tu)
	
	b=span/2.
	d=0.12*chord
	pegelamiet=calcAmietGershfeld(b,d,1.0,U,Tu*U,Lambda,freq,343.)
	pegelamiet2=calcAmiet(b,1.0,U,Tu*U,Lambda,freq,343.)
	
	print max(pegelamiet), Tu

	semilogx(freq,pegelamiet,'-',linewidth=1.0,color=collist[j],label=labellist[j])
	semilogx(freq,pegelamiet2,'--',linewidth=1.0,color=collist[j])

grid(color='0.5',linestyle=':',linewidth=0.2)
# savetxt('Tu_038_Amiet+Gershfeld.csv',np.column_stack((pegelamiet, pegelamiet2)), fmt='%5s', delimiter=',')
# savetxt('Tu_038_Amiet.csv',np.column_stack((freq, pegelamiet, pegelamiet2)), fmt='%5s', delimiter=',')
xticks((20,50,100,200,500,1000,2000,5000,10000,20000),('0,02','0,05','0,1','0,2','0,5','1','2','5','10','20'))
xlim(20,20000)
ylim(ymin,ymax)
legend(loc='lower left')
xlabel('$f_m$ in kHz',labelpad=1.0)
ylabel(r'$L_{p}$ in dB',labelpad=1.0)
gca().set_position([0.135,0.13,0.81,0.84] )
# savefig('Flatplate_LE_Tu038_.pdf',dpi=600)
savefig('naca0012_U5.pdf',dpi=600)
plt.show()