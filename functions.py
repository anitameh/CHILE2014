import numpy as np

import pyfits
import scipy
from scipy.stats import pearsonr


import math as ma
import matplotlib.pyplot as mp


import minuit


#----------------FUNCTIONS

imin=0
imax=0
jmin=0
jmax=0
data=0
B=0
Noise=0



def Source2D(x,y,xo,yo,sigmax,sigmay,rho,A):

	return A*ma.exp( (-(x-xo)**2.0/(2.0*sigmax**2.0)-(y-yo)**2.0/(2.0*sigmay**2.0)-rho*(x-xo)*(y-yo)/(sigmax*sigmay))/(1.0-rho**2.0)) /(2.0*ma.pi*sigmax*sigmay*ma.sqrt(1.0-rho**2.0))


def Source1D(x,xo,sigmax,A):

	return A*ma.exp(-(x-xo)**2.0/(2.0*sigmax**2.0))/ma.sqrt((2.0*ma.pi*sigmax))





def find(Data,Imin,Imax,Jmin,Jmax,background,noise,I):

	global imin
	imin=Imin
	global imax
	imax=Imax
	global jmin
	jmin=Jmin
	global jmax
	jmax=Jmax
	global data
	data=Data
	global B
	B=background
	global Noise
	Noise=noise	


	mi = minuit.Minuit(chi2,xo=(imin+imax)/2.0,yo=(jmin+jmax)/2.0,sigmax=1.5,sigmay=1.5,rho=0.0,A=I)
	mi.limits['sigmax'] = (1.0,50.0)
	mi.limits['sigmay'] = (1.0,50.0)
	mi.limits['rho']=(-1.0,1.0)
	#mi.values['rho']=0.0
	mi.limits['xo'] = (imin,imax)  
	mi.limits['yo'] = (jmin,jmax)  
	mi.limits['A']=(0.0,1e8)	
	#mi.values['A'] = I  
	#mi.printMode = 1
	results=np.zeros((5,2))
	try:
		mi.migrad()

		results=np.zeros(12)
		results[0],results[1]=mi.values['xo'],mi.errors['xo'] #x position
		results[2],results[3]=mi.values['yo'],mi.errors['yo'] #x position
		results[4],results[5]=mi.values['sigmax'],mi.errors['sigmax'] #Spread in x axis
		results[6],results[7]=mi.values['sigmay'],mi.errors['sigmay'] #Spread
		results[8],results[9]=mi.values['rho'],mi.errors['rho']

		results[10],results[11]=mi.values['A'],mi.errors['A'] # I

		return results
	except:
		return -np.ones(12)

def chi2(xo,yo,sigmax,sigmay,rho,A):
	
	chi=0.0
	
	for i in range(imin,imax):
		for j in range(jmin,jmax):
			chi=chi+(data[j,i]-Source2D(i,j,xo,yo,sigmax,sigmay,rho,A)-B)**2.0/(Noise**2.0)

	return chi



def find2(Data,Imin,Imax,Jmin,Jmax,background,noise):

	global imin
	imin=Imin
	global imax
	imax=Imax
	global jmin
	jmin=Jmin
	global jmax
	jmax=Jmax
	global data
	data=Data
	global B
	B=background
	global Noise
	Noise=noise

	mi = minuit.Minuit(chi22,x1o=(imax+imin)/2,y1o=(jmax+jmin)/2,x2o=(imax+imin)/2,y2o=(jmax+jmin)/2,sigma1x=1.5,sigma1y=1.5,sigma2x=1.5,sigma2y=1.5,A1=1e4,A2=1e4)


	mi.limits['sigma1x'] = (1.0,50.0)
	mi.limits['sigma1y'] = (1.0,50.0)
	mi.limits['sigma2x'] = (1.0,50.0)
	mi.limits['sigma2y'] = (1.0,50.0)
	mi.limits['x1o'] = (imin,imax)  
	mi.limits['y1o'] = (jmin,jmax)  
	mi.limits['x2o'] = (imin,imax)  
	mi.limits['y2o'] = (jmin,jmax) 
	mi.limits['A1'] =  (0,1e10) 
	mi.limits['A2'] = (0,1e10)  
	#mi.printMode = 1

	mi.migrad()

	results=np.zeros((10,2))
	results[0,:]=[mi.values['x1o'],mi.errors['x1o']] #x position
	results[1,:]=[mi.values['y1o'],mi.errors['y1o']] #x position
	results[2,:]=[mi.values['x2o'],mi.errors['x2o']] #x position
	results[3,:]=[mi.values['y2o'],mi.errors['y2o']] #x position

	results[4,:]=[mi.values['sigma1x'],mi.errors['sigma1x']] #Spread in x axis
	results[5,:]=[mi.values['sigma1y'],mi.errors['sigma1y']] #Spread in x axis
	results[6,:]=[mi.values['sigma2x'],mi.errors['sigma2x']] #Spread
	results[7,:]=[mi.values['sigma2y'],mi.errors['sigma2y']] #Spread
	results[8,:]=[mi.values['A1'],mi.errors['A1']] 
	results[9,:]=[mi.values['A2'],mi.errors['A2']] 
	return results



def chi22(x1o,y1o,x2o,y2o,sigma1x,sigma1y,sigma2x,sigma2y,A1,A2):
	
	chi=0.0

	for i in range(imin,imax):
		for j in range(jmin,jmax):
			chi=chi+(data[j,i]-Source2D(i,j,x1o,y1o,sigma1x,sigma1y,A1)-Source2D(i,j,x2o,y2o,sigma2x,sigma2y,A2)-B)**2.0/(Noise**2.0)

	return chi


def Goodness(D,model): #You give him just the block
	N=len(D[:,0])*len(D[0,:])
	x=np.zeros(N**2.0)
	y=np.zeros(N**2.0)
	for i in range(len(D[0,:])):
		for j in range(len(D[:,0])):

			x[i+N*j]=D[j,i]	
			y[i+N*j]=model[j,i]	
	return pearsonr(x,y)



def Goodness2(Data,jmin,jmax,imin,imax,xo,yo,sigmax,sigmay,rho,A,Back,noise):
	
	chi=0.0
	#Lx=len(Data[0,imin:imax])
	#Ly=len(Data[:,0])
	dif=np.zeros((jmax-jmin,imax-imin))
	for i in range(imax-imin):
		for j in range(jmax-jmin):
			dif[j,i]=data[jmin+j,imin+i]-Source2D(i+imin,j+jmin,xo,yo,sigmax,sigmay,rho,A)-Back
	mean=np.mean(dif)
	std=np.std(dif)

	for i in range(len(dif[0,:])):
		for j in range(len(dif[:,0])):
			chi=chi+(dif[j,i]-mean)**2.0/(std**2.0)

	return (chi/((jmax-jmin)*(imax-imin))-1.0)*ma.sqrt((jmax-jmin)*(imax-imin)/2.0)


def inter(i,j,q,ps1,ps2,Fin):

	f=0.0
	S=0.0
	a=2.0*ps2
	di=(i-N/2.0)*ps2
	dj=(j-N/2.0)*ps2
	
	ni=int(di/ps1+M/2.0-2.0*a/ps1)
	mi=int(dj/ps1+M/2.0-2.0*a/ps1)
	nmax=int(di/ps1+M/2.0+2.0*a/ps1)
	mmax=int(dj/ps1+M/2.0+2.0*a/ps1)
	if ni<0: ni=0
	if mi<0: mi=0
	if nmax<0: nmaxi=0
	if mmax<0: mmax=0
	if ni>M: ni=M-1
	if mi>M: mi=M-1
	if nmax>M: nmaxi=M-1
	if mmax>M: mmax=M-1
	for n in range(ni,nmax):

		dn=(n-M/2.0)*ps1
		k=0
		for m in range(mi,mmax+1):
			
			dm=(m-M/2.0)*ps1
			
			r=ma.sqrt((dn-di)**2.0+(dm-dj)**2.0)
			if r<a: 
				P=ma.cos(ma.pi*r/(2.0*a))
				f=f+P*Fin[q,n,m]*ps1**2.0
				S=S+P
				k=1
			elif k==1: break
	if S==0.0: return 0.0		
	else: return f/(S*ps2**2.0)

def smooth(data,a):
	Ni=len(data[0,:])
	Nj=len(data[:,0])

	F=np.zeros((Nj,Ni), dtype=np.float64)

	for i in range(Ni):

		for j in range(Nj):
				F[j,i]=interpol(i,j,a,data)
	return F

def interpol(i,j,a,data):
	S=0.0
	I=0.0
	kmin=max(int(i-a*2),0)
	kmax=min(int(i+a*2),len(data[0,:]))
	mmin=max(int(j-a*2),0)
	mmax=min(int(j+a*1),len(data[:,0]))
	for k in range(kmin,kmax):

		for m in range(mmin,mmax):

			r=ma.sqrt((m-j)**2.0+(k-i)**2.0)
			if r<a:
				P=ma.cos(ma.pi*r/(2.0*a))			
				I+=data[m,k]*P
				S+=P
	return I/S



def PixtoWC(i,j,CD1_1,CD1_2,CD2_1,CD2_2,CRPIX1,CRPIX2,PV,CRVAL1,CRVAL2):

	x=CD1_1*(i-CRPIX1)+CD1_2*(j-CRPIX2)
	y=CD2_1*(i-CRPIX1)+CD2_2*(j-CRPIX2)

	r=ma.sqrt(x**2.0+y**2.0)

	chi=PV[0,0]+PV[0,1]*x+PV[0,2]*y+PV[0,3]*r+PV[0,4]*x**2.0+PV[0,5]*x*y+PV[0,6]*y**2.0+PV[0,7]*x**3.0+PV[0,8]*x**2.0*y+PV[0,9]*x*y**2.0+PV[0,10]*y**3.0


	eta=PV[1,0]+PV[1,1]*y+PV[1,2]*x+PV[1,3]*r+PV[1,4]*y**2.0+PV[1,5]*x*y+PV[1,6]*x**2.0+PV[1,7]*y**3.0+PV[1,8]*y**2.0*x+PV[1,9]*y*x**2.0+PV[1,10]*x**3.0	


	alpha=ma.degrees( ma.atan( (ma.radians(chi)/ma.cos(ma.radians(CRVAL2)))/(1-ma.radians(eta)*ma.tan(ma.radians(CRVAL2)))) )
 
	RA=alpha+CRVAL1
	DEC= ma.degrees( ma.atan( (ma.radians(eta)+ma.tan(ma.radians(CRVAL2)))*ma.cos(ma.radians(alpha))/(1-ma.radians(eta)*ma.tan(ma.radians(CRVAL2) )) )  )

	return [RA,DEC]

# output: errPixtoWC calculates RA, DEC and their standard deviations
# input: delta_i, delta_j are very small values;
#        Sigma is a 2 by 2 matrix of the variance of (i,j), which is [[sigmax2,rho*sigmax*sigmay],[rho*sigmax*sigmay,sigmay2]]


def errPixtoWC(i,delta_i,j,delta_j,Sigma,CD1_1,CD1_2,CD2_1,CD2_2,CRPIX1,CRPIX2,PV,CRVAL1,CRVAL2):

        f_delta_x = PixtoWC(i+delta_i,j,CD1_1,CD1_2,CD2_1,CD2_2,CRPIX1,CRPIX2,PV,CRVAL1,CRVAL2)
        f_delta_y = PixtoWC(i,j+delta_j,CD1_1,CD1_2,CD2_1,CD2_2,CRPIX1,CRPIX2,PV,CRVAL1,CRVAL2)
        f = PixtoWC(i,j,CD1_1,CD1_2,CD2_1,CD2_2,CRPIX1,CRPIX2,PV,CRVAL1,CRVAL2)

        RA = f[0]
        DEC = f[1]

        RA_derivative = np.matrix( [[(f_delta_x[0] - RA)/delta_i, (f_delta_y[0] - RA)/delta_j]])
        DEC_derivative = np.matrix( [[(f_delta_x[1] - DEC)/delta_i, (f_delta_y[1] - DEC)/delta_j]])

        error_RA = np.sqrt(RA_derivative * Sigma * RA_derivative.T)
        error_DEC = np.sqrt(DEC_derivative * Sigma * DEC_derivative.T)
        
        return RA, DEC, error_RA, error_DEC



def isCosmicRay(sub_image,th, th_small):
	rows = len(sub_image)
	cols = len(sub_image[0])
	size_type = 'big'

	if(rows < 12 or cols < 12):
		th=th_small	
		size_type = 'small'

	if(rows > 25 or cols > 20):
		#take the center of the image
		i_center = int(rows /2)
		j_center = int(cols /2)
		sub_image = sub_image[i_center-10:i_center+11,j_center-10:j_center+11]

	rows = len(sub_image)
	cols = len(sub_image[0])

	sub_image_reshaped = np.reshape(sub_image,rows*cols,1)
	sub_image_kurtosis = scipy.stats.kurtosis(sub_image_reshaped,axis=0,fisher=False,bias=True)
	sub_image_skew     = scipy.stats.skew(sub_image_reshaped,axis=0,bias=True)

	if (sub_image_skew * sub_image_kurtosis) > th or sub_image_skew > 5 or sub_image_kurtosis > 20:
		return [ True , round(sub_image_kurtosis,3) , round(sub_image_skew,3) , round(sub_image_kurtosis*sub_image_skew,3) , size_type]
	else:
		return [ False, round(sub_image_kurtosis,3) , round(sub_image_skew,3) , round(sub_image_kurtosis*sub_image_skew,3) , size_type]




