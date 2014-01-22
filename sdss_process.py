import sys, copy, linecache
import numpy as np
import os
######################################################
#					FUNCTIONS
######################################################

def fun_distance(x1,y1,x2,y2,x1_err,y1_err,x2_err,y2_err,err_atro,threshold):
	err_x = np.sqrt( x1_err**2 + x2_err**2 + err_atro**2 )
	err_y = np.sqrt( y1_err**2 + y2_err**2 + err_atro**2 )
	distance = np.sqrt( ((x2-x1)**2 / err_x )  + (y2-y1)**2 / err_y )
	flag_threshold = False
	if (distance < threshold):
		flag_threshold = True
	return [distance, flag_threshold]
def fun_distanceLog(file1,file2):
	file1_count = 0
	file2_count = 0
	matches = []
	for row1 in file1.tolist():
		file2_count = 0
		for row2 in file2.tolist():
			x1 = row1[0]
			y1 = row1[1]
			x2 = row2[0]
			y2 = row2[2]
			x1_err = 0.0000
			y1_err = 0.0000
			x2_err = row2[1]
			y2_err = row2[3]
			astro_err = 0.000027
			threshold = 1.0
			dist = fun_distance(x1,y1,x2,y2,x1_err,y1_err,x2_err,y2_err,astro_err,threshold)
			if(dist[1]):
				row_result = np.zeros((1,4)).tolist()[0]
				row_result[0] = file1_count
				row_result[1] = file2_count
				row_result[2] = dist[0]
				row_result[3] = int(dist[1])
				matches.append(row_result)
			file2_count += 1
		file1_count += 1
	return matches
def fun_delete_tuple(list_tuples,obj1,obj2):
	result_list = []
	for row in list_tuples:
		if(row[0] !=  obj1 and row[1] !=  obj1) and (row[0] !=  obj2 and row[1] !=  obj2):
			result_list.append(list(row))
	return result_list
def fun_objCounter(catalog,size_f1,size_f2,ismatch):
	objects_counter = np.zeros((1 , size_f1+size_f2)).tolist()[0]
	for row in catalog:
		index1 = int(row[0])
		index2 = int(row[1])
		if(row[3] == ismatch):
			objects_counter[index1] += 1
			objects_counter[index2] += 1
	return objects_counter
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

######################################################
#					  MAIN
######################################################

#OPEN FILES
#>> python sdss_process.py /Users/dicotips/Desktop/Catalog_Calibration/sdss_file.csv /Users/dicotips/Desktop/Catalog_Calibration/input_file.csv /Users/dicotips/Desktop/Catalog_Calibration/output.cat

#file_sdss_path = str(sys.argv[1])
file_input_path = str(sys.argv[1])
#file_output_path = str(sys.argv[3])
file_output_path = '/home/local/school/student19/scratch/catalogue/SDSSmatched/' + os.path.basename(file_input_path)

file_sdss_path   = "/home/local/school/student19/scratch/src/sdss_file.csv"
#file_input_path  = "/Users/dicotips/Desktop/Catalog_Calibration/input_file.csv"
#file_output_path = "/Users/dicotips/Desktop/Catalog_Calibration/output.cat"

file_sdss = np.loadtxt(file_sdss_path,delimiter=',')
file_in   = np.loadtxt(file_input_path)

#Getting RA=, DEC=  from FILE1 and FILE2
file_sdss_coords = file_sdss[:,7:8+1]
file_in_coords = file_in[:,0:4]

#print file_in_coords


#Getting Matches
matches = fun_distanceLog(file_sdss_coords,file_in_coords)
result_matches_1_1  = fun_listMatches_1_1(matches, len(file_sdss),len(file_in))
result_matches_rank = fun_listMatches_Rank(matches, len(file_sdss),len(file_in))
result_matches = result_matches_1_1 + result_matches_rank

######################################################
#					  RESULTS
######################################################

#Index of columns U, U_Err, ObjId in SDSS_File
idx_u 		= int(9)
idx_u_err 	= int(14)

file_sdss_list = file_sdss.tolist()
file_in_list   = file_in.tolist()

#Generating Data for the OUTPUT_FILE
f = open(file_output_path, 'w')
for row in result_matches:
	u     = file_sdss_list[row[0]][idx_u]
	u_err = file_sdss_list[row[1]][idx_u_err]
	line_num_sdss = row[0]
	#extract N line data from files that matches
	line_f_in   = linecache.getline(file_input_path, row[1]+1)
	line_f_in   = line_f_in[0:len(line_f_in)-1].strip()
	line_f_sdss = linecache.getline(file_sdss_path , row[0]+1)
	line_f_sdss = line_f_sdss[0:len(line_f_in)-1].strip()
	#object_id_in   = line_f_in[0:line_f_in.index(',')].strip()
	object_id_sdss = line_f_sdss[0:line_f_sdss.index(',')].strip()
	object_id = object_id_sdss #[object_id_in,object_id_sdss]
	output_row = line_f_in + "\t" + str(u) + "\t" + str(u_err) + "\t" + str(line_num_sdss) +"\t"+ str(object_id) + "\n"
	f.write(output_row)
f.close()

