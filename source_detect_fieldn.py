import scipy, pyfits, copy, sys
import scipy.io as sci
import scipy.stats
import numpy as np
import os

import Queue
import functions as SbFun  #Sebastian Functions

#sys.setrecursionlimit(10000)



############################################################
############################################################

#Source Detection OUTPUT:  2D numpy.narray
def clean_data_func(data):
	rows = len(data)
	cols = len(data[0])

	data_size = [rows,cols]

	#Calculate column/row means and global means
	column_means = np.mean(data, axis=0)
	global_mean  = np.mean(column_means)

	#First remove the fake signals on the edge
	max_x = 1;
	max_y = 1;
	min_x = data_size[0];
	min_y = data_size[1];

	for i in range(0,data_size[0]):
		if(np.sum(data[i,:] > (global_mean - 5.0)) < 10):
			if(i < data_size[0] / 2):
				max_x = max(max_x,i)
			else:
				min_x = min(min_x,i)

	for i in range(0,data_size[1]):
		if(np.sum(data[:, i] > global_mean - 5.0) < 10):
			if(i < data_size[1] / 2):
				max_y = max(max_y, i)
			else:
				min_y = min(min_y, i)

	data_clean = data[(max_x+1):(min_x), (max_y+1):(min_y)];
	return data_clean

def detect_points_for_source(quantile_value_specified, data):
	data_size = [len(data),len(data[0])]
	bool_source_data_significant = copy.deepcopy(data)

	for i in range(0,data_size[0]):
		for j in range(0,data_size[1]):
			if (data[i,j] > quantile_value_specified):
				bool_source_data_significant[i, j] = 1;	
			else:
				bool_source_data_significant[i, j] = 0;
	return bool_source_data_significant

#Object Detection >>> OUTPUT:  [coords, 2D_list]

def propagate(image_org,i,j,objIndex):
	q_cache = Queue.Queue(maxsize=0)

	q_cache.put([i,j])
	while not q_cache.empty():
		actual_pixel = q_cache.get()
		p_i = actual_pixel[0]
		p_j = actual_pixel[1]

		size_i = len(image_org)
		size_j = len(image_org[0])

		if(image_org[p_i][p_j]==1):
			image_org[p_i][p_j] = objIndex
		elif(image_org[p_i][p_j]>1):
			continue

		if(p_i-1 < size_i and p_i-1 >= 0) and (p_j-1 < size_j and p_j-1 >= 0):
			if(image_org[p_i-1][p_j-1]==1):
				q_cache.put([p_i-1,p_j-1])	

		if(p_i-1 < size_i and p_i-1 >= 0) and (p_j < size_j and p_j >= 0):
			if(image_org[p_i-1][p_j]==1):
				q_cache.put([p_i-1,p_j])	

		if(p_i-1 < size_i and p_i-1 >= 0) and (p_j+1 < size_j and p_j+1 >= 0):
			if(image_org[p_i-1][p_j+1]==1):
				q_cache.put([p_i-1,p_j+1])

		if(p_i < size_i and p_i >= 0) and (p_j-1 < size_j and p_j-1 >= 0):
			if(image_org[p_i][p_j-1]==1):
				q_cache.put([p_i,p_j-1])
		if(p_i < size_i and p_i >= 0) and (p_j+1 < size_j and p_j+1 >= 0):
			if(image_org[p_i][p_j+1]==1):
				q_cache.put([p_i,p_j+1])

		if(p_i+1 < size_i and p_i+1 >= 0) and (p_j-1 < size_j and p_j-1 >= 0):
			if(image_org[p_i+1][p_j-1]==1):
				q_cache.put([p_i+1,p_j-1])	

		if(p_i+1 < size_i and p_i+1 >= 0) and (p_j < size_j and p_j >= 0):
			if(image_org[p_i+1][p_j]==1):
				q_cache.put([p_i+1,p_j])	

		if(p_i+1 < size_i and p_i+1 >= 0) and (p_j+1 < size_j and p_j+1 >= 0):
			if(image_org[p_i+1][p_j+1]==1):
				q_cache.put([p_i+1,p_j+1])

def findCoords(image, n):
	rows = len(image)
	cols = len(image[0])
	
	#Creating the result list of coords
	coords_row = [0]*6
	coords = []
	for i in range(0,n):
		coords.append(copy.deepcopy(coords_row))

	#Detection of the boxes coords [objIdx , x1, y1 , x2 , y2 , flag1, flag2]
	for i in range(0,rows):
		for j in range(0,cols):
			objIdx = image[i][j]

			if(objIdx <> 0):
				objIdx = image[i][j]
				coords[objIdx][0] = objIdx

				if(coords[objIdx][5] == 0):
					coords[objIdx][1] = i
					coords[objIdx][2] = j
					coords[objIdx][3] = i
					coords[objIdx][4] = j
					coords[objIdx][5] = 1

				#MIN_X
				if(i < coords[objIdx][1]):
					coords[objIdx][1] = i
				
				#MIN_Y				
				if(j < coords[objIdx][2]):
					coords[objIdx][2] = j
				
				#MAX_X
				if(i > coords[objIdx][3]):
					coords[objIdx][3] = i
				
				#MAX_Y
				if(j > coords[objIdx][4]):
					coords[objIdx][4] = j
	return coords
def printListM(mat,space):
	rows = len(mat)
	cols = len(mat[0])
	
	for i in range(0,rows):
		for j in range(0,cols):
			print str(mat[i][j])+space,
		print
def findObject(image_org):
	rows = len(image_org)
	cols = len(image_org[0])

	##CREATION OF A NEW IMAGE WITH
	image_result = copy.deepcopy(image_org)
	objIndex = 2

	for i in range(0,rows):
		for j in range(0,cols):
			if(image_result[i][j] == 1):
				propagate(image_result,i,j,objIndex)
				objIndex +=1

	##Standarazing number of objects list (by -1 for each object)
	for i in range(0,rows):
		for j in range(0,cols):
			if(image_result[i][j] > 0):
				image_result[i][j] = image_result[i][j]-1

	#Calculating Object Boxes coords
	coords = findCoords(image_result, objIndex-1)




	#printListM(image_org,'')
	#print
	#printListM(image_result,'')
	
	return [coords, np.array(image_result)]



#-------------------------MAIN-------------------------


# Load image data
file_path = sys.argv[1]
file_name = os.path.basename(file_path)

epoch = file_name.split('_')[-2]
field_num = file_name.split('_')[-4]
ccd_ns = file_name.split('_')[1]


fit= pyfits.open(file_path)
print 'processing %s' % file_path

# open output file

cat_out = '/home/local/school/student19/scratch/catalogue/field%s/image_%s_%s_%s.cat' % (field_num, field_num, ccd_ns, epoch)
fp_out = open(cat_out, 'a+') # Opens a file for both appending and reading.										Cat.append(Dat)


data 	= fit[0].data #[0:500,0:500] #extraer matriz de datos

header	= fit[0].header
t=header['EXPTIME']

data=data/t


CD1_1,CD1_2,CD2_1,CD2_2,CRPIX1,CRPIX2=header['CD1_1'],header['CD1_2'],header['CD2_1'],header['CD2_2'],header['CRPIX1'],header['CRPIX2']
PV=np.zeros((2,11))
PV[0,0]=header['PV1_0']
PV[0,1]=header['PV1_1']
PV[0,2]=header['PV1_2']
PV[0,3]=header['PV1_3']
PV[0,4]=header['PV1_4']
PV[0,5]=header['PV1_5']
PV[0,6]=header['PV1_6']
PV[0,7]=header['PV1_7']
PV[0,8]=header['PV1_8']
PV[0,9]=header['PV1_9']
PV[0,10]=header['PV1_10']

PV[1,0]=header['PV2_0']
PV[1,1]=header['PV2_1']
PV[1,2]=header['PV2_2']
PV[1,3]=header['PV2_3']
PV[1,4]=header['PV2_4']
PV[1,5]=header['PV2_5']
PV[1,6]=header['PV2_6']
PV[1,7]=header['PV2_7']
PV[1,8]=header['PV2_8']
PV[1,9]=header['PV2_9']
PV[1,10]=header['PV2_10']


CRVAL1=header['CRVAL1']
CRVAL2=header['CRVAL2']

data_cleaned = data

rows = len(data_cleaned)
cols = len(data_cleaned[0])
data_size = [rows,cols]

#Detect possible sources using a specified quantile

data_reshaped = np.reshape(data_cleaned,rows*cols,1);
quantiles_data = scipy.stats.mstats.mquantiles( data_reshaped, [0.9985, 0.998, 0.997, 0.993, 0.991, 0.987, 0.98, 0.96,0.8,0.7,0.6,0.5])

bool_source_data_significant_1 = detect_points_for_source(quantiles_data[0], data)



#Detect Objects from binary image
mat_list = bool_source_data_significant_1.tolist()
coords = findObject(mat_list)

subimage_4 = coords[1]

#printListM(coords[0],'\t')

#Calculating backgroun data
image1_list = data.tolist()
image3_list = bool_source_data_significant_1.tolist()

image_masked = np.ma.array(image1_list, mask=image3_list)

# plt.imshow(image3_list)
# plt.show()

global_mean	= np.mean(data)
global_std	= np.std(data)
background_mean	= image_masked.mean()
background_std	= image_masked.std()

print "gloal_mean:   " + str(global_mean)
print "gloab_stdev:  " + str(global_std)

print "backgr_mean:  " + str(background_mean)
print "backgr_stdev: " + str(background_std)


# print coords[0]
print len(coords[0])
Ny=len(data[:,0])
Nx=len(data[0,:])
print Nx, Ny
Cat=[]
stars=0

for i in range(1,len(coords[0])):


	y1 = coords[0][i][1]
	x1 = coords[0][i][2]
	y2 = coords[0][i][3]
	x2 = coords[0][i][4]

	if coords[0][i][4]-coords[0][i][2]>0 and coords[0][i][3]-coords[0][i][1]>0:

		imin=max(coords[0][i][2]-10,0)	
		imax=min(coords[0][i][4]+10,Nx)
		jmin=max(coords[0][i][1]-10,0)
		jmax=min(coords[0][i][3]+10,Ny)
		cos=SbFun.isCosmicRay(data[jmin:jmax,imin:imax],100,50)

		if cos[0]==False:	
			imin=max(coords[0][i][2]-5,0)	
			imax=min(coords[0][i][4]+5,Nx)
			jmin=max(coords[0][i][1]-5,0)
			jmax=min(coords[0][i][3]+5,Ny)

			I=np.sum(data[jmin:jmax,imin:imax]-background_mean) 
			r=SbFun.find(data,imin,imax,jmin,jmax,background_mean,background_std,I)


			if r[0]!=-1.0 and r[1]<5.0 and r[3]<5.0:

				model=np.zeros((jmax-jmin,imax-imin))

				g=SbFun.Goodness2(data,jmin,jmax,imin,imax,r[0],r[2],r[4],r[6],r[8],r[10],background_mean,background_std)
				np.append(r,g)

				if g<0.1:
					np.append(r,g)
					#Cat.append(r)
					stars+=1

					##asignar coordenadas
					sigma=np.zeros((2,2))
					sigma[0,0]=r[4]**2.0
					sigma[0,1]=r[8]*r[4]*r[6]

					sigma[1,0]=sigma[0,1]
					sigma[1,1]=r[6]**2.0

					b=(r[6]**2.0+r[4]**2.0)/2.0-np.sqrt(((r[6]**2.0-r[4]**2.0)/2.0)**2.0+(r[8]*r[6]*r[4])**2.0)

					a=(r[6]**2.0+r[4]**2.0)/2.0+np.sqrt(((r[6]**2.0-r[4]**2.0)/2.0)**2.0+(r[8]*r[6]*r[4])**2.0)

					r[0],r[2],r[1], r[3]= SbFun.errPixtoWC( r[0], r[1], r[2],r[3], sigma,CD1_1, CD1_2,CD2_1, CD2_2, CRPIX1, CRPIX2,PV,CRVAL1,CRVAL2)

					Dat=np.zeros(8)
					Dat[0:4]=r[0:4]
					Dat[4]=a
					Dat[5]=b
					Dat[6]=r[10]
					Dat[7]=r[11]
					Cat.append(Dat)
for i in xrange(len(Cat)):
	wrstr = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (Cat[i][0], Cat[i][1], Cat[i][2], Cat[i][3], Cat[i][4], Cat[i][5], Cat[i][6], Cat[i][7])
	fp_out.write(wrstr)



fp_out.close()



