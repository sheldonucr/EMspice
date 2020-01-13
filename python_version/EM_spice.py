import numpy as np
import os
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
	#os.system(cmd)	
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
	for i in range(1,(nx_total - 1)):
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
		#with open('Void.txt', 'w') as positionfile:
			#positionfile.writelines(' '.join(str(j) for j in i) + '\n' for i in Void_position) 
		#with open('Lvoid.txt', 'w') as lvfile:
			#lvfile.writelines(' '.join(str(j) for j in i) + '\n' for i in Void_length) 
					
	#clean_repo(tree_info, tree_width, tree_curden,tree_leng)
	stress_collection, Void_collection, Void_length_collection = simulation_process(tree_list,tree_stress,Void_position,Void_length,E,Z,cu_res,Omega)
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
def simulation_process(tree_list,tree_stress,Void_position,Void_length,E,Z,cu_res,Omega):
	
	stress_collection = []
	Void_collection = []
	Void_length_collection = []
		
	for i in range(len(tree_list)):
		new_stress, onetree_Void_position, onetree_Void_length = 	one_tree_simulation(tree_list[i],tree_stress[i],Void_position[i],Void_length[i],E,Z,cu_res,Omega,t_total)
		stress_collection.append(new_stress)
		Void_collection.append(onetree_Void_position)
		Void_length_collection.append(onetree_Void_length)
	#print(type(stress_collection[0][0]))
	
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
			onetree_Void_position = [left_branch_num] = 1
		else:
			onetree_Void_position[left_branch_num -1] = -1
			onetree_Void_position[left_branch_num] = 1
	return new_stress, onetree_Void_position, onetree_Void_length
def one_tree_simulation(one_tree,one_tree_stress,onetree_Void_position,onetree_Void_length,E,Z,cu_res,Omega,t_total):
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
		new_stress, onetree_Void_position, onetree_Void_length = solvegrotree(J_list_norm, one_tree,nx,t_total,kappa, Stress_norm,E,Z,cu_res,Omega,onetree_Void_position,onetree_Void_length)	
	
			
	return new_stress, onetree_Void_position, onetree_Void_length		
def solvegrotree(J_list_norm, one_tree,nx,t_total,kappa, Stress_norm,E,Z,cu_res,Omega,onetree_Void_position,onetree_Void_length):
	new_stress = []
	tstamp = normalize_time(t_total, one_tree,kappa)
	if onetree_Void_position[0] ==	1:
		print("t0")
		G,C,B = matrix_formation_growth(J_list_norm, one_tree,nx)
		pwl = [1] * len(B)

		sol_stress = backEuler_fdm(G,C,B,tstamp,Stress_norm,pwl)
		new_stress = denormalize_stress(one_tree, sol_stress,E,Z,cu_res,Omega)
			
	elif onetree_Void_position[-1] == -1:
		# inverse stress and tree_information
		Stress_norm_tmp, J_list_norm_tmp, tree_tmp = inverse_information(Stress_norm,J_list_norm,one_tree)
		G,C,B = matrix_formation_growth(J_list_norm_tmp, tree_tmp,nx)
		pwl = [1] * len(B)

		sol_stress_tmp = backEuler_fdm(G,C,B,tstamp,Stress_norm_tmp,pwl)
		new_stress_tmp = denormalize_stress(tree_tmp, sol_stress_tmp,E,Z,cu_res,Omega)
		new_stress , _, _ = inverse_information(new_stress_tmp ,J_list_norm_tmp,tree_tmp)
	
		print("t1")
	else:
		print("t2")
	return new_stress, onetree_Void_position, onetree_Void_length
def inverse_information(Stress_norm,J_list_norm,one_tree):
	Stress_norm_tmp = []
	J_list_norm_tmp = []
	for i in range(len(Stress_norm)):
		Stress_norm_tmp.append(Stress_norm[-1-i])
	for i in range(len(J_list_norm)):
		J_list_norm_tmp.append(J_list_norm[-1-i])
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
		branch_curden_tmp.append(branch_curden_origin[-1-i])

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
	G[0,0] = G[0,0] - 1
	G[0,1] = G[0,1] + 1
	G[-1,-1] = G[-1,-1] - 1
	G[-1,-2] = G[-1,-2] + 1
	G[0,0] = G[0,0] * (-1-(dx/delta))
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

def EM_simulation(T, D0, E, Ea, kB, Da, B0, Omega, cu_res, hCu, Ta_res, hTa, Z, kappa ,nx,wire_length, t_total,input_spice):
	onestep_simulation(T, D0, E, Ea, kB, Da, B0, Omega, cu_res, hCu, Ta_res, hTa, Z, kappa ,nx,wire_length, t_total,input_spice)	

if __name__ == '__main__':
	constant_file = "EM_spice_constant.txt"
	input_spice = "armcore.sp"	
	T, D0, E, Ea, kB, Da, B0, Omega, cu_res, hCu, Ta_res, hTa, Z, kappa ,nx,wire_length, t_total= get_constant(constant_file)
	EM_simulation(T, D0, E, Ea, kB, Da, B0, Omega, cu_res, hCu, Ta_res, hTa, Z, kappa ,nx,wire_length, t_total,input_spice )
		

