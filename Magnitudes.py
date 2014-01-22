import scipy, pyfits, copy, sys
import scipy.io as sci
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
import math as ma
import os

#input L,K, catalogues of a CCD

#output the same catalogues with magnitudes instead of counts and airmass

#open




for field_idx in xrange(40):
	print "processing field%s" % (field_idx+1)
	# input catalogues
	path1="/home/local/school/student19/scratch/catalogue/field%s/" % (field_idx+1)

	# output results
	path2="/home/local/school/student19/scratch/calibrated_results/field%s/" % (field_idx+1)
	
	if not os.path.exists(path2):
		os.makedirs(path2)

	Names_cat=os.listdir(path1)

	K=0.5
	eK=0.05

	Airmass=np.zeros(len(Names_cat))


	for i in xrange(len(Names_cat)):
		
		Cat=np.loadtxt(path1+Names_cat[i])
		# print i
	#-----------Assign Airmasses
		namedat=Names_cat[i].split("_")
		field=namedat[1]
		epoch=round(float(namedat[3][0:-4]), 5)

		pathair="/home/local/school/student19/scratch/src/airmass.list"
		arch=open(pathair, "r" )
		numlineas= len(arch.readlines( ))
		arch.seek(0)
		#print field, epoch

		for j in range(numlineas):

			linea=arch.readline()
			dat=str.split(linea)
		
			f=dat[1].split("_")
			ep=round(float(dat[2][0:-4]), 5)
	
			if ep==epoch and f[1] == str(field_idx + 1): #field 1 
				# print ep, epoch			
				# print "hello"
				a=float(dat[3])

				Airmass[i]=a
				break
		if Airmass[i]==0.0 and j==numlineas: 
			print i,Names_cat[i]
		arch.close()
		
	
	#-----------Assign Ls

		ccd=namedat[2]
		pathL="/home/local/school/student19/scratch/src/Ls.dat"
		arch=open(pathL, "r" )
		numlineas= len(arch.readlines( ))
		arch.seek(0)
		for j in range(numlineas):

			linea=arch.readline()
			dat=str.split(linea)


			if ccd==dat[0]:  

				L=float(dat[1])
				eL=float(dat[2])

				break
		arch.close()

		l=sum(1 for line in open(path1+Names_cat[i]))		
	
		for j in xrange(l):
			if l>1:
				C=Cat[j,6] #counts
				eC=Cat[j,7] #errcounts

				u=-2.5*ma.log10(C)+K*Airmass[i]+L
				eu=ma.sqrt((2.5*eC/(ma.log(10.0)*C))**2.0+eL**2.0+(eK*Airmass[i])**2.0)

				Cat[j,6]=u #magnitudes
				Cat[j,7]=eu #emagnitudes
			elif l==1:  
				C=Cat[6]
				eC=Cat[7]

				u=-2.5*ma.log10(C)+K*Airmass[i]+L
				eu=ma.sqrt((2.5*eC/(ma.log(10.0)*C))**2.0+eL**2.0+(eK*Airmass[i])**2.0)
				

		# write the results to disk
		f_handle = file(path2+Names_cat[i], 'a+')
		np.savetxt(f_handle,Cat)
		f_handle.close()

		
