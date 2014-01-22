import sys, io, copy
import numpy as np
import os
def fun_distance(x1,y1,x2,y2,x1_err,y1_err,x2_err,y2_err,err_atro,threshold):

	err_x = np.sqrt( x1_err**2 + x2_err**2 + err_atro**2 )
	err_y = np.sqrt( y1_err**2 + y2_err**2 + err_atro**2 )

	distance = np.sqrt( ((x2-x1)**2 / err_x**2 )  + (y2-y1)**2 / err_y**2 )

	flag_threshold = False
	if (distance < threshold):
		flag_threshold = True
	return [distance, flag_threshold]
def fun_distanceLog(file1,file2):
	file1_count = 0
	file2_count = 0

	match1 = []

	for row1 in file1.tolist():
		file2_count = 1
		for row2 in file2.tolist():
			x1 = row1[0]
			y1 = row1[2]
			x2 = row2[0]
			y2 = row2[2]
			x1_err = row1[1]
			y1_err = row1[3]
			x2_err = row2[1]
			y2_err = row2[3]
			astro_err = 0.000027
			threshold = 1.0 #----------------
		
			dist = fun_distance(x1,y1,x2,y2,x1_err,y1_err,x2_err,y2_err,astro_err,threshold)

			row_result = np.zeros((1,4)).tolist()[0]
			
			row_result[0] = file1_count
			row_result[1] = file2_count
			row_result[2] = dist[0]
			row_result[3] = int(dist[1])

			match1.append(row_result)

			file2_count += 1
		file1_count += 1

	return match1

def fun_objCounter(catalog,size_f1,size_f2,ismatch):
	objects_counter = np.zeros((1 , size_f1+size_f2)).tolist()[0]
	for row in catalog:
		index1 = int(row[0])
		index2 = int(row[1])
	
		if(row[3] == ismatch):
			objects_counter[index1] += 1
			objects_counter[index2] += 1

	return objects_counter
def fun_listMatches(catalog,size_f1, size_f2):

	objects_counter = fun_objCounter(catalog,size_f1,size_f2,1)

	matches_list = []

	#Select objects that match with >=1 object (using threshold)
	for objIndex in range(0,len(objects_counter)):
		#find object in catalog
		for row_catalog in catalog:
			if(row_catalog[0]==objIndex or row_catalog[1]==objIndex) and row_catalog[3]==1:
				matches_list.append(tuple(row_catalog))

		#Eliminate duplicates
		matches_list = list(set(matches_list))

		#Sort elements to display
		matches_list.sort()

	#Reducing and transforming the list [[source,destiny,distance],]
	matches_list_2=[]
	for row in matches_list:
		matches_list_2.append(list(row))

	#return catalog
	return matches_list_2
def fun_delete_tuple(list_tuples,obj1,obj2):
	result_list = []
	for row in list_tuples:
		if(row[0] !=  obj1 and row[1] !=  obj1) and (row[0] !=  obj2 and row[1] !=  obj2):
			result_list.append(list(row))
	return result_list

def fun_listMatches_1_1(matches_list,size_f1, size_f2):
	
	objects_counter_matches = fun_objCounter(matches_list,size_f1,size_f2,1)

	#Count left and right object connecions
	matches_list_2 = []

	for row in matches_list:
		if(objects_counter_matches[row[0]]==1 and objects_counter_matches[row[1]]==1):
			matches_list_2.append(row)

	return matches_list_2
def fun_listMatches_Rank(matches_list,size_f1, size_f2):
	
	objects_counter_matches = fun_objCounter(matches_list,size_f1,size_f2,1)

	#Count left and right object connecions
	matches_list_2 = []

	#Selecting matches n-n
	for row in matches_list:
		if(objects_counter_matches[row[0]]>1 or objects_counter_matches[row[1]]>1):
			matches_list_2.append(row)

	#Sorting elements by ASC distance. Convert to tuples each element
	matches_list_3 = []
	for row in matches_list_2:
		matches_list_3.append(tuple(row))

	matches_list_3.sort(key=lambda x: x[2])

	#Add as match  first element and objects related
	matches_list_near_k1 = []
	
	while matches_list_3:
		row = list(matches_list_3[0])
		matches_list_near_k1.append(row)
		matches_list_3 = fun_delete_tuple(matches_list_3, row[0],row[1])

	return matches_list_near_k1
def fun_list_no_matches(matches_1_1,matches_rank,size_f1,size_f2):
	matches_list = matches_1_1 + matches_rank
	no_matches_list = []

	#Generating Objects that not matched from Epoce1

	if len(matches_list)>0: var = np.array(matches_list)[:,0].tolist()
	else: var=[]
	for e1 in range(0,size_f1):
		if e1 not in var:
			no_matches_list.append([e1,None,None,0])

	#Generating Objects that not matched from Epoce2
	if len(matches_list)>0: 	var = np.array(matches_list)[:,1].tolist()
	else: var=[]

	for e2 in range(size_f1,size_f2+size_f2+1):
		if e2 not in var:
			no_matches_list.append([None,e2,None,0])
	return no_matches_list

def findMathes(file1,file2):
	#Calculate distance between objects in File1 and File2
	result_log = fun_distanceLog(file1, file2)
	result_matches = fun_listMatches(result_log, len(file1),len(file2))

	#Process matches 1-1
	result_matches_1_1 = fun_listMatches_1_1(result_matches, len(file1),len(file2))

	#Process matches N-N (top1-Rank minimizing distance)
	result_matches_rank = fun_listMatches_Rank(result_matches, len(file1),len(file2))

	#process no-mathces
	#result_no_matches = fun_list_no_matches(result_matches_1_1,result_matches_rank,len(file1),len(file2))

	#concatenating Lists of Matches and No-Matches
	total_list = result_matches_1_1 + result_matches_rank# + result_no_matches

	if len(total_list)>0: final_list = np.array(total_list)[:,0:2].tolist()
	else: final_list=[]
	return final_list

def fun_CatalogMatrix_isIn2(CatalogM, elem1, elem2):
	flag = -1
	index = 0
	for tkn in CatalogM:
		if(elem1 in tkn[1]) or (elem2 in tkn[1]):
			flag=index
			break
		index += 1
	return flag
def fun_CatalogMatrix_isIn1(CatalogM, elem1):
	flag = -1
	index = 0
	for tkn in CatalogM:
		if(elem1 in tkn[1]):
			flag=index
			break
		index += 1
	return flag
def fun_sumBoolArray(bool_array):
	suma = 0
	for i in bool_array:
		if (i):
			suma += 1
	return suma
def token_debug_repeated(token):
	file_num_list = map(lambda x: x[:x.index('_')],token)
	token_debugged = []
	dist=0
	index = 0
	for x in file_num_list:
		if(x not in file_num_list[:dist]):
			token_debugged.append(token[index])
		dist+=1
		index += 1
	#print token
	#print token_debugged
	return token_debugged

######################################################
#					MAIN
######################################################

field=sys.argv[1]
path1='/home/local/school/student19/scratch/calibrated_results/field'+field+"/" #input field i directory 
path2='/home/local/school/student19/scratch/lightcurves/field'+field +"/" #output field i directory

if not os.path.exists(path2):
	os.makedirs(path2)

print path1
files=os.listdir(path1)

#OPEN FILES
#args = sys.argv[1:-1]
#args_out = sys.argv[-1]


CCDs=['N1', 'N2', 'N3', 'N4', 'N5','N6','N7','N8', 'N9', 'N10', 'N11','N12','N13','N14','N15', 'N16','N17', 'N18', 'N19', 'N20','N21', 'N22','N23', 'N24','N25','N26', 'N27', 'N28', 'N29', 'N31','S1', 'S2', 'S3', 'S4', 'S5','S6','S7','S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 'S16','S17', 'S18', 'S19', 'S20','S21', 'S22','S23', 'S24', 'S25','S26', 'S27', 'S28', 'S29', 'S30', 'S31']

for ccd in CCDs:

	#Processing Matches interFiles
	args=[]	
	for f in files:
		if f.split('_')[2]==ccd: 
			args.append(f)
	print args
	if(len(args)>1):

		Catalog_Matrix = []
	
		flag_first_pass = False
		list_file = range(0,len(args))

		object_number = 0
		row = [0,[],[False]*len(args)]
		dist=1
		for idx_i in reversed(range(0,len(args))):
			for idx_j in range(0,idx_i):
				#print idx_j,
				#print idx_j+dist,
				#print flag_first_pass
				print "MATCHING:\tM_" + str(idx_j) + "_" + str(idx_j+dist) + "\t" + args[idx_j] + "\t" + args[idx_j+dist]

				file1 = np.loadtxt(path1+args[idx_j])
				file2 = np.loadtxt(path1+args[idx_j+dist])
				Matches = findMathes(file1,file2)

				counter = 1
				for i in Matches:
					if(i[0] is not None):
						obj1 = "%02d"%(idx_j)      + "_" + "%03d"%(int(i[0]))					
					else:
						obj1 = None
					
					if(i[1] is not None):
						obj2 = "%02d"%(idx_j+dist) + "_" + "%03d"%(int(i[1]))
					else:
						obj2 = None
				
					obj_idx = fun_CatalogMatrix_isIn2(Catalog_Matrix, obj1, obj2)
					if(obj_idx==-1):
						row1 = copy.deepcopy(row)
						row1[0] = object_number
						if(obj1 is not None):
							row1[1].append(obj1)
						if(obj2 is not None):
							row1[1].append(obj2)
						Catalog_Matrix.append(row1)
						object_number +=1
					else:
						if(obj1 is not None):
							if(obj1 not in Catalog_Matrix[obj_idx][1]):
								if(fun_CatalogMatrix_isIn1(Catalog_Matrix,obj1)==-1):
									Catalog_Matrix[obj_idx][1].append(obj1)
						if(obj2 is not None):
							if(obj2 not in Catalog_Matrix[obj_idx][1]):
								if(fun_CatalogMatrix_isIn1(Catalog_Matrix,obj2)==-1):
									Catalog_Matrix[obj_idx][1].append(obj2)
			flag_first_pass=True
			dist += 1

	#Processing the catalog to MATRIX putting true on them
	for row in Catalog_Matrix:
		for tkn in row[1]:
			file_idx = int(tkn[0:tkn.index('_')])
			if(row[2][file_idx] == False):
				row[2][file_idx] = True

	#Cleaning tokens
	for row in Catalog_Matrix:
		row[1]=token_debug_repeated(row[1])

	#Adding counter of object inter-files
	for row in Catalog_Matrix:
		row.append(fun_sumBoolArray(row[2]))

	#Discarding 1-time object in files
	obj_counter=0
	Catalog_Matrix_Final = []
	for row in Catalog_Matrix:
		if(row[-1] > 1):
			row1 = copy.deepcopy(row)
			row1[0] = obj_counter
			Catalog_Matrix_Final.append(row1)
			obj_counter += 1

	#Adding object number in files
	file_num = 0
	for _file in args:
		file_out_path = path2 +"/"+ ccd + "/" + _file

		f_in  = open(path1+_file,"r")

		if not os.path.exists(path2 +"/"+ ccd): os.makedirs(path2+"/"+ ccd)
		f_out = open(file_out_path,"w")
		f_size = len(f_in.readlines())
		f_in.seek(0)	
		row_num = 0
		object_in_file = []		#Objects included in file
		while row_num < f_size:
			coord = "%02d"%file_num + "_%03d"%row_num
			object_num = fun_CatalogMatrix_isIn1(Catalog_Matrix_Final,coord)
			row_file = f_in.readline()
			if(object_num != -1):
				row_final = str(object_num) + " " + _file.split('_')[3][:-4] + " " + row_file
				f_out.write(row_final)
				object_in_file.append(object_num)
			row_num += 1
		f_in.close()
		#Including NAN elements
		objects_total = range(0,len(Catalog_Matrix_Final))
		missing_objects = list(set(objects_total) - set(object_in_file))
		for obj in missing_objects:
			row_final = str(obj) + " nan nan nan nan nan nan nan nan nan \n"
			f_out.write(row_final)
		f_out.close()
		print "CREATED:" + _file
		file_num += 1


