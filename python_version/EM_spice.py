import numpy as np
import os
import matplotlib.pyplot as plt
from copy import deepcopy
class tree:
	def __init__(self,branch_num,branch_name, branch_length, branch_width,branch_curden):
		self.branch_num = branch_num
		self.branch_name = branch_name.copy()
		self.branch_length = branch_length.copy()
		self.branch_width = branch_width.copy()
		self.branch_curden = branch_curden.copy()
def get_constant(constant_file):
	with open (constant_file, 'r') as f:
		EM_constants = f.readlines()
		T = float(EM_constants[0].split(',')[-1])
		D0 = float(EM_constants[1].split(',')[-1])
		E = float(EM_constants[2].split(',')[-1])
		Ea_ratio = float(EM_constants[3].split(',')[-1])
		Ea = Ea_ratio * E
		kB = float(EM_constants[4].split(',')[-1])
		Da = D0 * np.exp(-Ea/(kB * T))
		B0 = float(EM_constants[5].split(',')[-1])
		Omega = float(EM_constants[6].split(',')[-1])
		cu_res = float(EM_constants[7].split(',')[-1])
		hCu = float(EM_constants[8].split(',')[-1])
		Ta_res = float(EM_constants[9].split(',')[-1])
		hTa = float(EM_constants[10].split(',')[-1])
		Z = float(EM_constants[11].split(',')[-1])
		kappa = (Da * B0 * Omega) / (kB*T)
		nx = float(EM_constants[12].split(',')[-1]) # the number of pieces of fdm
		wire_length = float(EM_constants[13].split(',')[-1])
		#t_total = float(EM_constants[14].split(',')[-1])# normalized total time
		#tstep = float(EM_constants[15].split(',')[-1]) # normalized tstep
		#tstamp = [] # every time step point
		t_total = float(EM_constants[14].split(',')[-1]) # total time for a big step

		
	return T, D0, E, Ea, kB, Da, B0, Omega, cu_res, hCu, Ta_res, hTa, Z, kappa,nx,wire_length, t_total  
def run_emcmd(input_spice):
	cmd = "./em_cmd " + input_spice
	os.system(cmd)	
def get_all_width(tree_width):
	all_branch_width =[] 
	with open(tree_width, 'r') as f:
		all_branch_width_origin = f.readlines()
		for branch_width_origin in all_branch_width_origin:
			branch_width = branch_width_origin[1:]
			all_branch_width.append(float(branch_width))
	return all_branch_width 
def get_all_curden(tree_curden):
	all_branch_curden = []
	with open(tree_curden, 'r') as f:
		all_branch_curden_origin = f.readlines()
		for branch_curden_origin in all_branch_curden_origin:
			branch_curden = branch_curden_origin.split()[-1]
			all_branch_curden.append(float(branch_curden))
	return all_branch_curden 
def get_all_leng(tree_leng):
	all_branch_leng = []
	with open(tree_leng, 'r') as f:
		all_branch_leng_origin = f.readlines()
		for branch_leng_origin in all_branch_leng_origin:
			branch_leng = branch_leng_origin[1:]
			all_branch_leng.append(float(branch_leng))
	return all_branch_leng 	 
def gettree(tree_info,tree_width, tree_curden, tree_leng):
	all_branch_width = get_all_width(tree_width)
	all_branch_curden = get_all_curden(tree_curden)
	all_branch_leng = get_all_leng(tree_leng)
	tree_list = gettree_info(tree_info, all_branch_width, all_branch_curden, all_branch_leng)
	print(tree_list[0].branch_width )	
	return tree_list
def get_one_branch_name(line):
	name = line.split()[0]
	return name
def gettree_info(tree_info, all_branch_width, all_branch_curden, all_branch_leng):
	tree_list  = []
	branch_name = [] # branch name of a tree
	branch_num = 0 # number of branch of a tree
	branch_counter = 0 # branch counter of the file		
	branch_length = [] 
	branch_width = []
	branch_curden = []
	with open (tree_info, 'r') as f:
		lines = f.readlines()
		for line in lines:
			if line.split()[0] == 'xx':
				newtree = tree(branch_num,branch_name, branch_length, branch_width,branch_curden)
				tree_list.append(newtree)
				branch_name = [] # branch name of a tree
				branch_num = 0 # number of branch of a tree
				branch_length = [] 
				branch_width = []
				branch_curden = []
			else:	
				one_branch_name = get_one_branch_name(line)
				branch_name.append(one_branch_name)
				branch_num = branch_num + 1
				branch_length.append(all_branch_leng[branch_counter])
				branch_width.append(all_branch_width[branch_counter])
				branch_curden.append(all_branch_curden[branch_counter])
				branch_counter = branch_counter + 1
	return 	tree_list			
#def get_tree_width(tree_width,tree_num):
#	brchwidth = [0] * tree_nnum
#	brch
#def init_parameters(treenum,initC,void,Lvoid, Lvoid1,Lvoid2,max_stress,max_stress_location,nx_total):
#	for i in range (treenum):
#		initC.append([])
#		void.append(0)
#		Lvoid.append(0)
#		Lvoid1.append(0)
#		Lvoid2.append(0)
#		max_stress.append(0)
#		max_stress_location.append(0)
#		nx_total.append(nx * brchnum[i]) - (nrchnum(i) - 1)
#		for j in range(nx_total[i]):
#			initC[i].append(0)
#		initC.append(nx_total)
#	Lvoid_last = Lvoid.copy
#	Lvoid_last1 = Lvoid1.copy
#	Lvoid_last2 = Lvoid2.copy
#
#	return initC, void, Lvoid, Lvoid1, Lvid2, max_stress, max_stress_location, nx_total, initC, Lvoid_last, Lvoid_last1, Lvoid_last2
#def simulation_process(init_flag, treenum, initC,void,Lvoid, Lvoid1,Lvoid2,max_stress,max_stress_location,nx_total, iterIdx):
#	if init_flag == 1:
#		initC = []
#		void = []
#		Lvoid = []
#		Lvoid1 = []
#		Lvoid2 = []
#		max_stress = []
#		max_stress_location = []
#		nx_total = []
#		init_parameters(treenum,initC,void,Lvoid, Lvoid1,Lvoid2,max_stress,max_stress_location,nx_total)
#		init_flag = 0
	# file name to store stress, curden and Lvoid
#	fname1 = 'u_stress_' + str(iterIdx) + '.txt'
#	fname2 = 'u_curden_' + str(iterIdx) + '.txt'
#	fname3 = 'u_Lvoid_' + str(iterIdx) + '.txt'
#	fid_resis = 'u_resis.txt'
#	for i in range(treenum):
#		Lvoid0 = 0
#		Lvoid01 = 0
#		Lvoid02 = 0
#		Lvoid_temp = 0
		#for k in range(brchnum[i])
#		FEM_run = 1
#		for j in range(FEM_run):
#			if void[i] == 0:
#				G,C,B = matrix_formation(nx,wires[i], J_list_norm[i])	
def matrix_formation_nuc(J_list_norm0, one_tree,nx):
	
	J_list_norm = np.reshape(np.asarray(J_list_norm0),(len(J_list_norm0),1))
	
	n_wires = one_tree.branch_num
		
	w_lengs = np.asarray(one_tree.branch_length)
	num_Junctions = n_wires - 1
	nx_total = int((nx * n_wires) - num_Junctions)
	
	G = np.zeros([nx_total,nx_total])
	B = np.zeros([nx_total,n_wires])
	for i in range(1,(nx_total - 2)):
		G[i,i-1] = G[i,i-1] + 1
		G[i,i] = G[i,i] - 2
		G[i,i + 1] = G[i,i + 1] + 1
	G[0,0] = G[0,0] - 1
	G[0,1] = G[0,1] + 1
	G[-1,-1] = G[-1,-1] - 1
	G[-1,-2] = G[-1,-2] + 1
	index = 0
	t = np.zeros([1,n_wires])
	for i in range(n_wires):
		B[index,i] = 1
		B[int(index + nx - 1), i] = -1
		index = int(index + nx - 1)
	C = np.ones([nx_total])
	dx = 1/(nx_total - 1)
	tmp =np.divide(J_list_norm,dx)
	B = np.dot(B,tmp)
	G = -G
	G = np.divide(G,(dx*dx))
	return G, C, B
	
def onestep_simulation(T, D0, E, Ea, kB, Da, B0, Omega, cu_res, hCu, Ta_res, hTa, Z, kappa ,nx,wire_length, t_total,input_spice):
	run_emcmd(input_spice)
	tree_info = "tree_info.txt"
	tree_width = "tree_width.txt"
	tree_curden = "tree_curden.txt"
	tree_leng = "tree_leng.txt"
	stress_file = 'u_stress.txt'
	Void_file = 'Void.txt'
	Lvoid_file = 'Lvoid.txt'	
	tree_list = gettree(tree_info, tree_width, tree_curden, tree_leng) # tree_list is a list of tree structure.
	tree_stress = []
	if os.path.isfile(stress_file ):
		print ("have initial stress")
		u_stress = stress_file 	
		tree_stress = get_stress_from_file(u_stress)
			
	else:
		tree_stress = getinitialstreee(tree_list, nx)
	if os.path.isfile(Void_file):
		Void_position, Void_length = getinitialvoid(Void_file ,Lvoid_file )
	else:
		Void_position = [] # for each tree, number of zero = branch number 0 means no void, 1 means void is on the left, -1 means on the write
		Void_length = [] # length of void on the branch, if 0 means no void.	
		for i in range (len(tree_list)):
			Void_position.append([])
			Void_length.append([])
			for j in range (tree_list[i].branch_num):
				Void_position[i].append(0)
				Void_length[i].append(0)
					
	clean_repo(tree_info, tree_width, tree_curden,tree_leng)
	Void_length_ini = deepcopy(Void_length)
	stress_collection, Void_collection, Void_length_collection = simulation_process(tree_list,tree_stress,Void_position,Void_length,E,Z,cu_res,Omega,B0)
	update_sp(Void_collection,Void_length_collection,tree_list,input_spice,Void_length_ini,Ta_res,hTa,hCu)	
	update_file(stress_collection, Void_collection, Void_length_collection, stress_file, Void_file,Lvoid_file,tree_list,nx)
def update_file(stress_collection, Void_collection, Void_length_collection, stress_file, Void_file,Lvoid_file,tree_list,nx):
	with open(stress_file, 'w') as stressfile:
		for i in range(len(tree_list)):
			for j in range(len(stress_collection[i])):
				if j%(nx-1) == 0 and j!=len(stress_collection[i]) - 1:
					b_num = int(j/(nx-1))
					b_name = tree_list[i].branch_name[b_num]
					stress_data = sum(stress_collection[i][j].tolist())
					outputline = b_name + ' ' + str(stress_data) + '\n'
					stressfile.writelines(outputline)
				else:
					stress_data = sum(stress_collection[i][j].tolist())
	
					outputline = ' ' + str(stress_data) + '\n'
					stressfile.writelines(outputline)

		#stressfile.writelines(' '.join(str(sum(j.tolist())) for j in i) + '\n' for i in stress_collection) 
	with open(Void_file, 'w') as positionfile:
		positionfile.writelines(' '.join(str(j) for j in i) + '\n' for i in Void_collection) 
	with open(Lvoid_file, 'w') as lvfile:
		lvfile.writelines(' '.join(str(j) for j in i) + '\n' for i in Void_length_collection) 

def update_sp(Void_position_collection,Void_length_collection,tree_list,input_spice,Void_length_origin,Ta_res,hTa,hCu):
	failed_branch, failed_void_length, failed_crit_length,origin_void_length_b, failed_position = get_growth_tree(Void_position_collection,Void_length_collection,tree_list,Void_length_origin)# get trees whose void is larger than cross section
	
	print("FC")
	print(failed_crit_length)	
	early_failed_branch,early_failed_void_length,early_failed_crit_length,late_failed_branch,late_failed_void_length,late_failed_crit_length,late_failed_origin_length,early_position, late_position =  early_late_failure(failed_branch, failed_void_length, failed_crit_length,origin_void_length_b,failed_position)
	modify_sp(early_failed_branch, early_failed_void_length, early_failed_crit_length, late_failed_branch, late_failed_void_length, late_failed_crit_length,input_spice,late_failed_origin_length,Ta_res,hTa,hCu,early_position, late_position)	


def early_late_failure(failed_branch, failed_void_length, failed_crit_length,Void_length_origin, failed_position):
	void_list = []
	early_failed_branch = [] 
	early_failed_void_length = [] 
	early_failed_crit_length = [] 
	late_failed_branch = [] 
	late_failed_void_length = [] 
	late_failed_crit_length = []
	late_failed_origin_length = []
	early_position = [] 
	late_position = []
	for i in range(len(failed_branch)):
		# R 25 is earlyt failure and others are late. 
		if failed_branch[i][0:3] == "R25": 
			early_failed_branch.append(failed_branch[i]) 
			early_failed_void_length.append(failed_void_length[i]) 
			early_failed_crit_length.append(failed_crit_length[i])
			early_position.append(failed_position[i])
		else:
			late_failed_branch.append(failed_branch[i]) 
			late_failed_void_length.append(failed_void_length[i]) 
			late_failed_crit_length.append(failed_crit_length[i])
			late_failed_origin_length.append(Void_length_origin[i])
			late_position.append(failed_position[i]) 

	return 	early_failed_branch,early_failed_void_length,early_failed_crit_length,late_failed_branch,late_failed_void_length,late_failed_crit_length,late_failed_origin_length,early_position, late_position

def modify_sp(early_failed_branch, early_failed_void_length, early_failed_crit_length, late_failed_branch, late_failed_void_length, late_failed_crit_length,input_spice,late_failed_origin_length,Ta_res,hTa,hCu,early_position, late_position):
	with open (input_spice, 'r') as f:
		all_line = f.readlines()
		modified_list = [] # list of modified branches with branch id
		for i in range (len(all_line)):	
			if all_line[i].split()[0] in early_failed_branch:
				cur_line = all_line[i]
				cur_line_sp2 =  cur_line.split()
				line1_id = i
				line2_id = i + 1
				if cur_line_sp2[0] in modified_list:
					continue
				index = early_failed_branch.index(all_line[i].split()[0])
				# early failure need to modify two lines.
				next_line = ""
				ref_line = ""
				# find if it is on the begining or end
				if  early_position[index] == 1:
			 	
					pre_line1 = all_line[i -1]
					pre_line_sp1 = pre_line1.split()
					if pre_line_sp1[1] != cur_line_sp2[2] and pre_line_sp1[2] != cur_line_sp2[1]:
						ref_line =  all_line[i + 1]
						modified_list.append(cur_line_sp2[0] )
					else:
						
						next_line = cur_line
						cur_line =  pre_line1
						line1_id = i -1
						line2_id = i
						modified_list.append(cur_line_sp2[0] )
						modified_list.append(pre_line_sp1 [0] )	
				else:
					pre_line1 = all_line[i + 1]
					pre_line_sp1 = pre_line1.split()				
					if pre_line_sp1[1] != cur_line_sp2[2] and pre_line_sp1[2] != cur_line_sp2[1]:
						ref_line =  all_line[i - 1]
						modified_list.append(cur_line_sp2[0] )
					else:
						next_line =  pre_line1
						modified_list.append(cur_line_sp2[0] )
						modified_list.append(pre_line_sp1 [0] )	
	
				early_failure_line1,early_failure_line2, failed_point_location = get_early_failure_line(early_failed_branch[index], early_failed_void_length[index], early_failed_crit_length[index],cur_line,early_position[index],next_line,ref_line)
				writing_location = "v(" + failed_point_location + ")"
				all_line[-2] = all_line[-2].strip() + writing_location + "\n"
				if early_failure_line2 == "":
					print(line1_id)
					print(early_failure_line1)
					print(all_line[line1_id])
					all_line[line1_id] = early_failure_line1
				else:
					print(line1_id)
					print(early_failure_line1)
					print(all_line[line1_id])
					print(line2_id)
					print(early_failure_line2)
					print(all_line[line2_id])	
					all_line[line1_id] = early_failure_line1 + "\n"
					all_line[line2_id] = early_failure_line2 + "\n"
	
		
			if all_line[i].split()[0] in late_failed_branch:
				index = late_failed_branch.index(all_line[i].split()[0])
				late_failure_line = get_late_failure_line(late_failed_branch[index], late_failed_void_length[index], late_failed_crit_length[index],all_line[i],late_failed_origin_length[index],Ta_res,hTa,hCu,late_position[index])
				print(i)
				print(late_failure_line)
				print(all_line[i])
				all_line[i] = late_failure_line + "\n"
	write_new_sp(all_line, input_spice)
def write_new_sp(all_line, input_spice):
	os.remove(input_spice)
	with open(input_spice, 'w') as newsp:
		for new_line in all_line:
			newsp.write(new_line)

def get_early_failure_line(early_branch, early_void_length, early_crit_length,origin_line,position,next_line,ref_line):
	newline1 = ""
	newline2 = ""
	
	origin_line_sp = origin_line.split()
	next_line_sp = next_line.split()
	if next_line == "":
		ref_position1 = ref_line.split()[1]
		ref_position2 = ref_line.split()[2]
		
		
		old_id_l = origin_line_sp[1]
		x_coor_l = int(old_id_l.split("_")[1])
		y_coor_l = int(old_id_l.split("_")[2])
		old_id_r = origin_line_sp[2]
		x_coor_r = int(old_id_r.split("_")[1])
		y_coor_r = int(old_id_r.split("_")[2])
		if x_coor_l ==  x_coor_r:
			if old_id_l == ref_position1 or old_id_l == ref_position2:
				new_id_coor = y_coor_r + 1
				new_part = old_id_l.split("_")[0] + "_" + old_id_l.split("_")[1] + "_" + str(new_id_coor)
				newline1 = origin_line_sp[0] + " " + origin_line_sp[1] + " " + new_part + " " + origin_line_sp[3]
				
				newline2 = ""
			elif old_id_r == ref_position1 or old_id_r == ref_position2:
				new_id_coor = y_coor_l + 1
				new_part = old_id_l.split("_")[0] + "_" + old_id_l.split("_")[1] + "_" + str(new_id_coor)
				newline1 = origin_line_sp[0] + " " + new_part + " " + origin_line_sp[2] + " " + origin_line_sp[3]
				newline2 = ""
		else:
			if old_id_l == ref_position1 or old_id_l == ref_position2:
				new_id_coor = x_coor_r + 1
				new_part = old_id_l.split("_")[0] + "_" + str(new_id_coor) + "_" + old_id_l.split("_")[2]
				newline1 = origin_line_sp[0] + " " + origin_line_sp[1] + " " + new_part + " " + origin_line_sp[3]
			elif old_id_r == ref_position1 or old_id_r == ref_position2:
				new_id_coor = x_coor_l + 1
				new_part = old_id_l.split("_")[0] + "_" + str(new_id_coor) + "_" + old_id_l.split("_")[1]
				newline1 = origin_line_sp[0] + " " + new_part + " " + origin_line_sp[2] + " " + origin_line_sp[3]
				newline2 = ""
		
	else:
		origin_line1 = origin_line_sp[1]
		origin_line2 = origin_line_sp[2]

		next_line1 = next_line_sp[1]
		next_line2 = next_line_sp[2]
		x_coor_l = int(origin_line1.split("_")[1])
		y_coor_l = int(origin_line1.split("_")[2])
		x_coor_r = int(origin_line2.split("_")[1])
		y_coor_r = int(origin_line2.split("_")[2])
		nx_coor_l = int(next_line1.split("_")[1])
		ny_coor_l = int(next_line1.split("_")[2])
		nx_coor_r = int(next_line2.split("_")[1])
		ny_coor_r = int(next_line2.split("_")[2])
		if x_coor_l ==  x_coor_r:
			if origin_line1 == next_line1:
				new_id_coor = y_coor_l + 1
				new_part = origin_line1 .split("_")[0] + "_" + origin_line1.split("_")[1] + "_" + str(new_id_coor)
				newline1 = origin_line_sp[0] + " " + new_part+ " " + next_line_sp[2] + " " + origin_line_sp[3]
				newline2 = next_line_sp[0] + " " +  new_part+ " " + next_line_sp[2] + " " + next_line_sp[3]
			elif origin_line1 == next_line2:
				new_id_coor = y_coor_l + 1
				new_part = origin_line1 .split("_")[0] + "_" + origin_line1.split("_")[1] + "_" + str(new_id_coor)
				newline1 = origin_line_sp[0] + " " + new_part+ " " + origin_line_sp[2] + " " + origin_line_sp[3]
				newline2 = next_line_sp[0] + " " + next_line_sp[1] + " " + new_part + " " + next_line_sp[3]
			elif origin_line2 == next_line1:
				new_id_coor = y_coor_r + 1

				new_part = origin_line1 .split("_")[0] + "_" + str(new_id_coor) + "_" + origin_line1.split("_")[2] 
				newline1 = origin_line_sp[0] + " " + next_line_sp[1] + " " + new_part + " " + origin_line_sp[3]
				newline2 = next_line_sp[0] + " " +  new_part+ " " + next_line_sp[2] + " " + next_line_sp[3]
			elif origin_line2 == next_line2:
				new_id_coor = y_coor_r + 1

				new_part = origin_line1 .split("_")[0] + "_" + str(new_id_coor) + "_" + origin_line1.split("_")[2]
				newline1 = origin_line_sp[0] + " " + next_line_sp[1] + " " + new_part + " " + origin_line_sp[3] 
				newline2 = next_line_sp[0] + " " + next_line_sp[1] + " " + new_part + " " + next_line_sp[3]
			else:
				pass
		
		
	return newline1, newline2, new_part
def get_deltR(late_void_length,late_failed_origin_length, failed_branch_width,Ta_res,hTa,hCu):
	print("deltR data")
	print(late_void_length )
	print(late_failed_origin_length)
	print(failed_branch_width)
	delt_R = Ta_res*(late_void_length - late_failed_origin_length)/(hTa*failed_branch_width)
	print(delt_R)
	return delt_R
def get_late_failure_line(late_branch, late_void_length, late_crit_length,origin_line,late_failed_origin_length,Ta_res,hTa,hCu,position):
	# here critical length is width
	delt_R = get_deltR(late_void_length,late_failed_origin_length,late_crit_length,Ta_res,hTa,hCu)
	origin_line_sp = origin_line.split()
	new_line = origin_line_sp[0] + " "+ origin_line_sp[1] + " " + origin_line_sp[2] + " " + str(float(origin_line_sp[3]) + delt_R)
	return new_line 
def get_growth_tree(Void_position,Void_length,tree_list,Void_length_origin):
	failed_branch = []
	failed_void_length = []
	failed_crit_length = []
	origin_void_length_b = [] # origin void_length on branch
	origin_void_length_f = [] # origin void length of failed branch
	all_position = []
	failed_position = []
	for i in range(len(tree_list)):
	
		if all(v == 0 for v in Void_position[i]):
			pass
		else:
			crit_void = []
			void_length = []
			branch_id = []	
			if Void_position[i][0] == 1:
				crit_void.append(tree_list[i].branch_width[0])
				void_length.append(Void_length[i][0])
				branch_id.append(tree_list[i].branch_name[0])
				origin_void_length_b.append(Void_length_origin[i][0])
				all_position.append(1) 
			elif Void_position[i][-1] == -1:
				crit_void.append(tree_list[i].branch_width[-1])
				void_length.append(Void_length[i][-1])
				branch_id.append(tree_list[i].branch_name[-1])
				origin_void_length_b.append(Void_length_origin[i][-1])
				all_position.append(-1)  
			else:
				left_branch = Void_position[i].index(-1)
				right_branch = Void_position[i].index(1)
				print("RL name")
				print(tree_list[i].branch_name[left_branch])
				print(tree_list[i].branch_name[right_branch])
				crit_void.append(tree_list[i].branch_width[left_branch])
				crit_void.append(tree_list[i].branch_width[right_branch])
				
				void_length.append(Void_length[i][left_branch ])
				void_length.append(Void_length[i][right_branch ])
				branch_id.append(tree_list[i].branch_name[left_branch ])	
				branch_id.append(tree_list[i].branch_name[right_branch ])
				origin_void_length_b.append(Void_length_origin[i][left_branch ])
				origin_void_length_b.append(Void_length_origin[i][right_branch]) 
				all_position.append(-1)
				all_position.append(1)  
			print("CV") 
			print(crit_void)
			print(void_length)
			for i in range(len(crit_void)):
				if void_length[i] > crit_void[i]:
					failed_branch.append(branch_id[i])
					failed_void_length.append(void_length[i])
					failed_crit_length.append(crit_void[i])	
					origin_void_length_f.append(origin_void_length_b[i])			
					failed_position.append(all_position[i])	
	return 	failed_branch, failed_void_length, failed_crit_length, origin_void_length_f, failed_position
def simulation_process(tree_list,tree_stress,Void_position,Void_length,E,Z,cu_res,Omega,B0):
	
	stress_collection = []
	Void_collection = []
	Void_length_collection = []
			
	for i in range(len(tree_list)):
		new_stress, onetree_Void_position, onetree_Void_length = 	one_tree_simulation(tree_list[i],tree_stress[i],Void_position[i],Void_length[i],E,Z,cu_res,Omega,t_total,B0)
		stress_collection.append(new_stress)
		Void_collection.append(onetree_Void_position)
		Void_length_collection.append(onetree_Void_length)
	
	return stress_collection, Void_collection, Void_length_collection			 
	#J_list_norm = normalize_vurden(tree_list) # all the current density over its maximum absolute value

def solvenuctree(J_list_norm, one_tree,nx,t_total,kappa, Stress_norm,E,Z,cu_res,Omega,onetree_Void_position,onetree_Void_length):
	G,C,B = matrix_formation_nuc(J_list_norm, one_tree,nx)
	tstamp = normalize_time(t_total, one_tree,kappa)
	pwl = [1] * len(B)
	sol_stress = backEuler_fdm(G,C,B,tstamp,Stress_norm,pwl)
	new_stress = denormalize_stress(one_tree, sol_stress,E,Z,cu_res,Omega)
	max_new_stress = max(new_stress)
	if max_new_stress > 5e8:
		position = new_stress.index(max(new_stress))
		left_branch_num = int(position//(nx - 1))
		if (left_branch_num) == len(onetree_Void_position):
			onetree_Void_position[left_branch_num -1] = -1
		elif left_branch_num == 0:
			onetree_Void_position [left_branch_num] = 1
		else:
			onetree_Void_position[left_branch_num -1] = -1
			onetree_Void_position[left_branch_num] = 1
	
	#plt.plot(new_stress)
	#plt.show()
	return new_stress, onetree_Void_position, onetree_Void_length
def one_tree_simulation(one_tree,one_tree_stress,onetree_Void_position,onetree_Void_length,E,Z,cu_res,Omega,t_total,B0):
	J_list_norm = normalize_curden(one_tree)
	Stress_norm = normalize_stress(one_tree, one_tree_stress,E,Z,cu_res,Omega)
	if all(v == 0 for v in onetree_Void_position):
		new_stress, onetree_Void_position, onetree_Void_length = solvenuctree(J_list_norm, one_tree,nx,t_total,kappa, Stress_norm,E,Z,cu_res,Omega,onetree_Void_position,onetree_Void_length)	
		#G,C,B = matrix_formation_nuc(J_list_norm, one_tree,nx)
		#tstamp = noramlize_time(t_total, one_tree,kappa)
		#pwl = [1] * len(B)
	
		#sol_stress = backEuler_fdm(G,C,B,tstamp,Stress_norm,pwl)
		#new_stress = denormalize_stress(one_tree, sol_stress,E,Z,cu_res,Omega)
		#max_new_stress = max(new_stress)
		#if max_new_stress > 5e8:
		#	position = new_stress.index(max(new_stress))
		#	left_branch_num = int(position//(nx - 1)) 
		#	if (left_branch_num) == len(onetree_Void_position) - 1:
		#		onetree_Void_position[left_branch_num] = -1
		#	else:
		#		onetree_Void_position[left_branch_num] = -1
		#		onetree_Void_position[left_branch_num + 1] = 1


	else:
		new_stress, onetree_Void_position, onetree_Void_length = solvegrotree(J_list_norm, one_tree,nx,t_total,kappa, Stress_norm,E,Z,cu_res,Omega,onetree_Void_position,onetree_Void_length,B0)	
	
			
	return new_stress, onetree_Void_position, onetree_Void_length		
def solvegrotree(J_list_norm, one_tree,nx,t_total,kappa, Stress_norm,E,Z,cu_res,Omega,onetree_Void_position,onetree_Void_length,B0):
	new_stress = []
	tstamp = normalize_time(t_total, one_tree,kappa)
	if onetree_Void_position[0] ==	1:
		G,C,B = matrix_formation_growth(J_list_norm, one_tree,nx)
		pwl = [1] * len(B)

		sol_stress = backEuler_fdm(G,C,B,tstamp,Stress_norm,pwl)
		new_stress = denormalize_stress(one_tree, sol_stress,E,Z,cu_res,Omega)
		L_Void  = sum(calculate_void(new_stress,one_tree,B0).tolist())

		onetree_Void_length[0] = L_Void	
	elif onetree_Void_position[-1] == -1:
		# inverse stress and tree_information
		Stress_norm_tmp, J_list_norm_tmp, tree_tmp = inverse_information(Stress_norm,J_list_norm,one_tree)
		G,C,B = matrix_formation_growth(J_list_norm_tmp, tree_tmp,nx)
		pwl = [1] * len(B)

		
		sol_stress_tmp = backEuler_fdm(G,C,B,tstamp,Stress_norm_tmp,pwl)
		new_stress_tmp = denormalize_stress(tree_tmp, sol_stress_tmp,E,Z,cu_res,Omega)
		new_stress , _, _ = inverse_information(new_stress_tmp ,J_list_norm_tmp,tree_tmp)
		
		
		L_Void  = sum(calculate_void(new_stress,one_tree,B0).tolist())
	
		onetree_Void_length[-1] = L_Void
	else:
		#print(sum(Stress_norm))
		# cut the wire on the void
		Stress_norm_1, J_list_norm_1, tree_1, Stress_norm_2, J_list_norm_2, tree_2  =  cuttree(onetree_Void_position, Stress_norm,J_list_norm,one_tree,nx)
		
		G1,C1,B1 = matrix_formation_growth(J_list_norm_1, tree_1,nx)
		pwl_1 = [1] * len(B1)
		#print(B1)
		#plt.plot(Stress_norm_1)
		#plt.show()	
		#print("case2" + str(Stress_norm_1[0]))
		sol_stress_1 = backEuler_fdm(G1,C1,B1,tstamp,Stress_norm_1,pwl_1)
		#plt.plot(sol_stress_1)
		#plt.show()
		#print(sum(sol_stress_1))
		G2,C2,B2 = matrix_formation_growth(J_list_norm_2, tree_2,nx)
		pwl_2 = [1] * len(B2)
		#plt.plot(Stress_norm_2)
		#plt.show()	
		#print(Stress_norm_2[0])
		sol_stress_2 = backEuler_fdm(G2,C2,B2,tstamp,Stress_norm_2,pwl_2)
		#plt.plot(sol_stress_2)
		#plt.show()
		#print(sum(sol_stress_2))
		com_solstress = combine_stress(sol_stress_1, sol_stress_2)
		new_stress = denormalize_stress(one_tree, com_solstress,E,Z,cu_res,Omega)
		#plt.plot(new_stress)
		#plt.show()
		Lvoid_left, Lvoid_right = cal_void_midcase(tree_1,tree_2,new_stress,B0,nx)
		if Lvoid_left < 0:
			Lvoid_left = 0
		if Lvoid_right < 0:
			Lvoid_right = 0
		for i in range(len(onetree_Void_position)):
			if onetree_Void_position[i] == -1:
				onetree_Void_length[i] = Lvoid_left
			if onetree_Void_position[i] == 1:
				onetree_Void_length[i] = Lvoid_right
			
	return new_stress, onetree_Void_position, onetree_Void_length

def cal_void_midcase(tree_1,tree_2,new_stress,B0,nx):
	left_branch_num =  tree_1.branch_num
	right_branch_num = tree_2.branch_num
	left_node_num = int(left_branch_num * (nx - 1) + 1)
	right_node_num = int(right_branch_num * (nx - 1) + 1)
	stress_1 = new_stress[0:left_node_num]
	stress_2 = new_stress[left_node_num - 1 : len(new_stress)]
	Lvoid_left = sum(calculate_void(stress_1,tree_1,B0).tolist())
	Lvoid_right = sum(calculate_void(stress_2,tree_2,B0).tolist())
	#print(Lvoid_left)
	#print(Lvoid_right)
	return Lvoid_left, Lvoid_right 	
def calculate_void(new_stress,one_tree,B0):
	n_wires = one_tree.branch_num
	L_total = sum(one_tree.branch_length)
	num_Junctions = n_wires - 1
	nx_total = int((nx * n_wires) - num_Junctions)
	dx = L_total/(nx_total - 1)
	Stress_integration = 0
	for i in range(len(new_stress)-1):
		Stress_integration = Stress_integration - (new_stress[i] + new_stress[i + 1]) /2 * dx
	L_Void = Stress_integration / B0
	return L_Void	
def combine_stress(sol_stress_1, sol_stress_2):
	com_solstress = []
	for i in range (len(sol_stress_1)):
		com_solstress.append(sol_stress_1[-1-i])
	com_solstress[-1] = (sol_stress_1[0] + sol_stress_2[0])/2
	for i in range (len(sol_stress_2) - 1):
		com_solstress.append(sol_stress_2[1 + i])
	return com_solstress	
def cuttree(onetree_Void_position, Stress_norm,J_list_norm,one_tree,nx):
	left_branch = onetree_Void_position.index(1)
	void_position = left_branch * (nx -1)
	Stress_norm_1_tmp = []
	Stress_norm_1 = []
	Stress_norm_2 = []
	J_list_norm_1_tmp = []
	J_list_norm_1 = []
	J_list_norm_2 = []
	for i in range(int(void_position + 1)):
		Stress_norm_1_tmp.append(Stress_norm[i])
	for i in range (int(void_position), len(Stress_norm),1):
		Stress_norm_2.append(Stress_norm[i])
	for i in range	(int(left_branch)):
		J_list_norm_1_tmp.append(J_list_norm[i])
	for i in range(int(left_branch), len(J_list_norm),1):
		J_list_norm_2.append(J_list_norm[i])
	branch_num_tmp = one_tree.branch_num
	branch_name_origin = one_tree.branch_name
	branch_length_origin = one_tree.branch_length
	branch_width_origin = one_tree.branch_width
	branch_curden_origin = one_tree.branch_curden
	branch_name_1 = []
	branch_length_1 = []
	branch_width_1 = []
	branch_curden_1 = []
	branch_name_2 = []
	branch_length_2 = []
	branch_width_2 = []
	branch_curden_2 = []

	branch_num_1 = left_branch
	branch_num_2 = branch_num_tmp - branch_num_1	
	for i in range (branch_num_1):
		branch_name_1.append(branch_name_origin[i])
		branch_length_1.append(branch_length_origin[i])
		branch_width_1.append(branch_width_origin[i])
		branch_curden_1.append(branch_curden_origin[i])
	for i in range (branch_num_1,branch_num_tmp,1):
		branch_name_2.append(branch_name_origin[i])
		branch_length_2.append(branch_length_origin[i])
		branch_width_2.append(branch_width_origin[i])
		branch_curden_2.append(branch_curden_origin[i])

	tree_1_tmp = tree(branch_num_1,branch_name_1, branch_length_1, branch_width_1,branch_curden_1)

	tree_2 = tree(branch_num_2,branch_name_2, branch_length_2, branch_width_2,branch_curden_2)
	Stress_norm_1, J_list_norm_1, tree_1 = inverse_information(Stress_norm_1_tmp,J_list_norm_1_tmp,tree_1_tmp)
	return Stress_norm_1, J_list_norm_1, tree_1, Stress_norm_2, J_list_norm_2, tree_2 				
def inverse_information(Stress_norm,J_list_norm,one_tree):
	# note that J has positive and negative value, which represent the direction of current, so the sign of j also need to be flipped. 
	Stress_norm_tmp = []
	J_list_norm_tmp = []
	for i in range(len(Stress_norm)):
		Stress_norm_tmp.append(Stress_norm[-1-i])
	for i in range(len(J_list_norm)):
		J_list_norm_tmp.append(-J_list_norm[-1-i])
	branch_num_tmp = one_tree.branch_num
	branch_name_origin = one_tree.branch_name
	branch_length_origin = one_tree.branch_length
	branch_width_origin = one_tree.branch_width
	branch_curden_origin = one_tree.branch_curden
	branch_name_tmp = []
	branch_length_tmp = []
	branch_width_tmp = []
	branch_curden_tmp = []
	for i in range(len(branch_name_origin)):
		branch_name_tmp.append(branch_name_origin[-1-i])
		branch_length_tmp.append(branch_length_origin[-1-i])
		branch_width_tmp.append(branch_width_origin[-1-i])
		branch_curden_tmp.append(-branch_curden_origin[-1-i])

	tree_tmp = tree(branch_num_tmp,branch_name_tmp, branch_length_tmp, branch_width_tmp,branch_curden_tmp)
	return 	Stress_norm_tmp, J_list_norm_tmp, tree_tmp 

def matrix_formation_growth(J_list_norm0, one_tree,nx):
	delta = 3e-7
	J_list_norm = np.reshape(np.asarray(J_list_norm0),(len(J_list_norm0),1))
	
	n_wires = one_tree.branch_num
	
	w_lengs = np.asarray(one_tree.branch_length)
	num_Junctions = n_wires - 1
	nx_total = int((nx * n_wires) - num_Junctions)
	dx = 1/(nx_total - 1)
	G = np.zeros([nx_total,nx_total])
	B = np.zeros([nx_total,n_wires])
	for i in range(1,(nx_total - 1)):
		G[i,i-1] = G[i,i-1] + 1
		G[i,i] = G[i,i] - 2
		G[i,i + 1] = G[i,i + 1] + 1
	#G[0,0] = G[0,0] - 1
	G[0,1] = G[0,1] + 1
	G[-1,-1] = G[-1,-1] - 1
	G[-1,-2] = G[-1,-2] + 1
	G[0,0] = (-1-(dx/delta))
	index = 0
	t = np.zeros([1,n_wires])
	for i in range(n_wires):
		B[index,i] = 1
		B[int(index + nx - 1), i] = -1
		index = int(index + nx - 1)
	C = np.ones([nx_total])
	tmp =np.divide(J_list_norm,dx)
	B = np.dot(B,tmp)
	B [0] = 0
	G = -G
	G = np.divide(G,(dx*dx))
	return G, C, B
def normalize_time(t_total, one_tree,kappa):
	onetree_branchleng = one_tree.branch_length.copy()
	treeleng = sum(onetree_branchleng)
	nor_t_total = t_total /(treeleng*treeleng) * kappa
	tstep = nor_t_total / 10
	tstamp = []
	for i in range(int(nor_t_total/tstep)):
		tstamp.append(i*tstep)
	return tstamp
def denormalize_stress(one_tree, sol_stress,E,Z,cu_res,Omega):
	onetree_curden = one_tree.branch_curden.copy()
	max_curden = abs(max(onetree_curden,key = abs)) 
	g1 = (E*Z*cu_res*max_curden)/Omega  
	onetree_branchleng = one_tree.branch_length.copy()
	treeleng = sum(onetree_branchleng)
	denormalized_stress = []
	for i in range(len(sol_stress)):
		denormalized_stress.append(sol_stress[i]*g1*treeleng)
	return denormalized_stress	
def backEuler_fdm(G,C,B,tstamp,Stress_norm,pwl):
	C_fdm = np.diag(C)
	G_fdm = np.asarray(G)
	B_fdm = np.asarray(B)
	Stress_norm_fdm = np.asarray(Stress_norm)
	Stress_norm_fdm = np.reshape(Stress_norm_fdm, (len(Stress_norm_fdm),1))
	u_fdm = np.reshape(np.asarray(pwl),(len(pwl),1))

	dt = tstamp[2] - tstamp[1]
	for i in range (len(tstamp) - 1):
		item1 = np.divide(C_fdm,dt)
		
		item2 = np.add(G_fdm, item1)
		item3 = np.multiply(B,u_fdm[i+1])
		item4 = np.dot(item1,Stress_norm_fdm)
		item5 = np.add(item3,item4)
		item6 = np.linalg.lstsq(item2,item5)[0]
		Stress_norm_fdm = np.copy(item6) 	
		#Stress_norm_fdm = (G_fdm + C_fdm / dt)\ (B * u_fdm[i + 1] + (C / dt) * Stress_norm_fdm)
	return Stress_norm_fdm
def normalize_stress(one_tree,one_tree_stress,E,Z,cu_res,Omega):
	onetree_curden = one_tree.branch_curden.copy()
	max_curden = abs(max(onetree_curden,key = abs)) 
	g1 = (E*Z*cu_res*max_curden)/Omega  
	onetree_branchleng = one_tree.branch_length.copy()
	treeleng = sum(onetree_branchleng)
	normalized_stress = []
	for i in range(len(one_tree_stress)):
		normalized_stress.append(one_tree_stress[i]/g1/treeleng)
	return normalized_stress	 
def normalize_curden(one_tree):
	onetree_curden = one_tree.branch_curden.copy()
	J_list_norm  = []
	max_curden = abs(max(onetree_curden,key = abs)) 
	for i in range(len(onetree_curden)):
		J_list_norm.append(onetree_curden[i]/max_curden)	
	return J_list_norm
def getinitialvoid(voidfile, Lvoidfile):
	Void_position = []
	Void_length = []
	with open (voidfile, 'r') as f:
		line = f.readline()
		while line:
			one_tree_position = []
			line_sp = line.split()
			for num in line_sp:
				one_tree_position.append(int(num))
			Void_position.append(one_tree_position)
			line = f.readline()
	with open (Lvoidfile, 'r') as f:
		line = f.readline()
		while line:
			one_tree_length = []
			line_sp = line.split()
			for num in line_sp:
				one_tree_length.append(float(num))
			Void_length.append(one_tree_length)
			line = f.readline()	
	return Void_position, Void_length 
def get_stress_from_file(u_stress):
	tree_stress = []
	with open (u_stress, 'r') as f:
		line_num = 0
		line = f.readline()
		line_sp = line.split()
		name = line_sp[0]
		tree_name = name.split("-")[0] + "-" + name.split("-")[1]
		
		current_tree_name = tree_name
		one_tree_stress = []
		c = 0
		while line:
			line_sp = line.split()
			if len(line_sp) == 2:
				name = line_sp[0]
				tree_name = name.split("-")[0] + "-" + name.split("-")[1]
				if tree_name != current_tree_name:
					c = c + 1
					tree_stress.append(one_tree_stress)
					one_tree_stress = []
					current_tree_name = tree_name
				stress = float(line_sp[1])
				one_tree_stress.append(stress)
		 
				
			else:
				stress = float(line_sp[0])
				
				one_tree_stress.append(stress)
			line = f.readline()	
		
		tree_stress.append(one_tree_stress)
	return tree_stress				
def getinitialstreee(tree_list,nx):
	# set all 0 for initial stress
	tree_stress = []
	for i in range(len(tree_list)):
		one_tree_stress = []
		branch_num = tree_list[i].branch_num
		stress_num = int(branch_num * (nx -1) + 1)
		for i in range(stress_num):
			one_tree_stress.append(0) # add initial stress for nodes.
		tree_stress.append(one_tree_stress)	
	return tree_stress	
def clean_repo(tree_info, tree_width, tree_curden,tree_leng):
	os.remove(tree_info)
	os.remove(tree_width)
	os.remove(tree_curden)
	os.remove(tree_leng)

def EM_simulation(T, D0, E, Ea, kB, Da, B0, Omega, cu_res, hCu, Ta_res, hTa, Z, kappa ,nx,wire_length, t_total,input_spice,simulation_number):
	for i in range(simulation_number):
		onestep_simulation(T, D0, E, Ea, kB, Da, B0, Omega, cu_res, hCu, Ta_res, hTa, Z, kappa ,nx,wire_length, t_total,input_spice)	

if __name__ == '__main__':
	constant_file = "EM_spice_constant.txt"
	input_spice = "armcore.sp"
	simulation_number = 1	
	T, D0, E, Ea, kB, Da, B0, Omega, cu_res, hCu, Ta_res, hTa, Z, kappa ,nx,wire_length, t_total= get_constant(constant_file)
	EM_simulation(T, D0, E, Ea, kB, Da, B0, Omega, cu_res, hCu, Ta_res, hTa, Z, kappa ,nx,wire_length, t_total,input_spice, simulation_number)
		

