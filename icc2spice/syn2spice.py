import numpy as np
import csv as cv
#import function as fun 

class branch:

	def __init__(self,resid, res_value,v_drop, current, cur_den, layer, width, start, end,wiretype):
		self.resid = resid
		self.res_value = res_value
		self.v_drop = v_drop
		self.current = current
		self.cur_den = cur_den
		self.layer = layer
		self.width = width
		self.start = start
		self.end = end
		self.wiretype = wiretype

	def get_data(self,resid, res_value,v_drop, current, cur_den, layer, width, start, end):
		self.resid = resid
		self.res_value = res_value
		self.v_drop = v_drop
		self.current = current
		self.cur_den = cur_den
		self.layer = layer
		self.width = width
		self.start = start
		self.end = end
	def get_type(self,wiretype):
		self.wiretype = wiretype
def read_line(layout_name):
	all_layer = []
	layer_info = []
	raw_layout = open(layout_name,encoding='utf-8')
	branch_info = []
	for line in raw_layout:

		if(line[0:12]=='* layer_data'):
			line_list=line.split()
			all_layer.append(int(line_list[2]))
			layer_info.append(line_list[4])
			
		if(line[0]=='#'):
			line_list=line.split( ) 
       #new_line is a list=[ id, res_value, drop_value(mV), current_value(mA), em_value (current_density)(A/cm), layer, width, x, y, x1, y1, ] length unit is nm
			new_line=[int(line_list[1]), float(line_list[2])/1000,float(line_list[5]),float(line_list[6])/1000,float(line_list[7]),int(line_list[8]),float(line_list[9]), int(line_list[10]),int(line_list[11]),int(line_list[12]),int(line_list[13])]
			branch_info.append(new_line)
	return branch_info, all_layer,layer_info

def reform_branch(raw_branch,all_layer, layer_info):
	branch_list = []
	for i in range(len(raw_branch)):
		id = raw_branch[i][0]
		res_value = raw_branch[i][1]
		v_drop = raw_branch[i][2]
		current = raw_branch[i][3]
		cur_den = raw_branch[i][4]
		layer = raw_branch[i][5]
		width = raw_branch[i][6]
		start = (raw_branch[i][7], raw_branch[i][8])
		end= (raw_branch[i][9], raw_branch[i][10])
		type_index = all_layer.index(layer)
		wiretype = layer_info[type_index][0]
		#new_branch = branch()
		new_branch = branch(id, res_value,v_drop, current, cur_den, layer, width, start, end,wiretype)
		#new_branch.get_data(id, res_value,v_drop, current, cur_den, layer, width, start, end)
	
		branch_list.append(new_branch)
	return(branch_list)

def seperate_layer(branch_list):
	layer_name = []
	branch_layer = []
	layer_type = []
	for branch in branch_list:
		if branch.layer in layer_name:
			branch_index = layer_name.index(branch.layer)
			branch_layer[branch_index].append(branch)
			
		
		else:
			layer_name.append(branch.layer)
			layer_type.append(branch.wiretype)
			branch_layer.append([branch])
	return layer_name, branch_layer,layer_type
def tree_position(current_bl, direc):
	tree = []
	cord_list = []
	sorted_branch = []
	for branch in current_bl:
		
		if direc == 0:
			direc_1 = 1
		elif direc == 1:
			direc_1 = 0
		cord_list.append(branch.start[direc_1])
	branch_order = np.argsort(cord_list)		
	for i in branch_order:
		sorted_branch.append(current_bl[i])
	for branch in sorted_branch:
		if len(tree) == 0:
			tree.append(branch)
		else:
			last_branch = tree[-1]
			if last_branch.end[direc_1] == branch.start[direc_1]:
				tree.append(branch)
			else:
				tree.append(branch)				 
				
	return tree

def get_currenttree(current_bl):
	current_tree = []
	layer_b = [] # seperate branches to different ordinate
	direc = -1# direc 0 for 0, direc 1 for 1

	index_list = [] # 0 for x coordinate and 1 for y coordinate
	if (current_bl[0].start[0] == current_bl[0].end[0]):
		direc = 0
		
		for branch in current_bl:
			if branch.start[0] not in index_list:
				index_list.append(branch.start[0])
				layer_b.append([])

				
				
			
	elif (current_bl[0].start[1] == current_bl[0].end[1]):	
		direc= 1
		for branch in current_bl:
			if branch.start[1] not in index_list:
				index_list.append(branch.start[1])
				layer_b.append([])
	
	for branch in current_bl:
		
		lb_index = index_list.index(branch.start[direc])
		layer_b[lb_index].append(branch)
	
	for current_bl in layer_b:
		one_tree = tree_position(current_bl, direc)
		current_tree.append(one_tree)
		
	return current_tree

	
def get_tree(layer_name, branch_layer,layer_type):
	tree_list = [];
	tree_layername = []
	for i in range(len(layer_name)):
		if layer_type[i] == 'v':
			continue
		else:
			
			current_layer = layer_name[i]
			layer_b = branch_layer[i]
			layer_tree = get_currenttree(layer_b)	
			tree_list.append(layer_tree)
			
	return tree_list, tree_layername


def collect_nodes(branch_list):
	node_position = []
	node_current = [] 
	for branch in branch_list:
		start_current  = branch.current
		end_current = 0 - branch.current 
		if branch.start not in node_position:
			node_position.append(branch.start)
			node_current.append(start_current)
		else:
			node_index_start = node_position.index(branch.start)
			node_current[node_index_start] + start_current
		if branch.end not in node_position:
			node_position.append(branch.end)
			node_current.append(end_current)
		else:
			node_index_end = node_position.index(branch.end)
			node_current[node_index_end] + end_current
	
		
	return node_position, node_current
def write_voltagesource(voltage_tree): #layer 41, connect with pad
	#node is named by their position for example (10,20 is named by N_10_20)
	pad_node = []
	for tree in voltage_tree:
		i = 0
		changed = 0
		
		while changed == 0 and i < len(tree):
			if tree[i].current > 0:
				pad_node.append(tree[0].start)
				changed = 1	
			elif tree[i].current < 0:
				pad_node.append(tree[-1].end)
				changed = 1
			else: 
				i = i + 1	
	# format Voltage source name, node name, 0, vdd(1.8)
	# node is named by its position like n_x_y
	voltage_source = []
	for i in range (len(pad_node)):
		voltagename = 'v' + str(i)
		
		nodename = 'n_' + str(pad_node[i][0]) + '_' + str(pad_node[i][1])
		vss = 0
		vdd = 1.8
		voltage_source.append([voltagename, nodename, str(vss), str(vdd)])	
	#write trees 
			
	return voltage_source
#	with open(filename, "a") as outputfile:
		
#def write_powergrid():
#def write_currentsource(): 

def write_tree(tree_layer,tree_layer_name):
	tree_output = []
	via_id = [i for i, x in enumerate(layer_type) if x == "v"]
	print(len(tree_layer))
	for i in range(len(tree_layer)):
		layer = tree_layer[i]# one layer
		layer_id = tree_layer_name[i]
		
		for j in range(len(tree_layer[i])): #one tree 
			for k in range (len(tree_layer[i][j])):
				branch = tree_layer[i][j][k]
				branch_name = 'R' + str(layer_id) + '_'+ str(j) + '_' + str(k) 
				node_name_start = 'n' + str(layer_id) + '_' +str(branch.start[0]) + '_' +str(branch.start[1]) 	
				node_name_end = 'n' + str(layer_id) + '_' + str(branch.end[0])+ '_' + str(branch.end[1])
				#node_name_start = 'n' + '_' +str(branch.start[0]) + '_' +str(branch.start[1]) 	
				#node_name_end = 'n' + '_' + str(branch.end[0])+ '_' + str(branch.end[1])

				res_name = str(branch.res_value)		
				tree_output.append([branch_name,  node_name_start,  node_name_end , res_name])
	return tree_output

def write_Csource(nodes_position,node_current):
	current_output = []
	for i in range(len(nodes_position)):
		cur_nmae = 'iB00_'+ str(i) + '_v'
		node_name = 'n1_'+str(nodes_position[i][0]) + '_'+ str(nodes_position[i][1])
		node_cur_name = str(node_current[i])
		current_output.append([cur_nmae, node_name , '0' ,  node_cur_name])
	return current_output

def write_spice(layer_name, branch_layer, tree_layer, layer_type, nodes_position, node_current,filename ): 
	# layer 41 is MRDL, connect with voltage source
	voltage_layerid = layer_name.index(41)
	voltage_layer = branch_layer[voltage_layerid]
	via_id = [i for i, x in enumerate(layer_type) if x == "v"]
	tree_id = voltage_layerid
	for inx in via_id:
		if inx< voltage_layerid:
			tree_id = tree_id - 1
	voltage_tree = tree_layer[tree_id]
	
	voltage_source_output = write_voltagesource(voltage_tree)
	tree_layer_name = []
	for i in range(len(layer_name)):
		if layer_type[i] == 'm':
			tree_layer_name.append(layer_name[i])
		
	tree_output = write_tree(tree_layer, tree_layer_name)
	current_output = write_Csource(nodes_position,node_current)
	output_list=voltage_source_output +  tree_output + current_output + [['.end']]
	with open(filename, "w") as f:
		writer = cv.writer(f)
		writer.writerows(output_list) 
	with open(filename, "r") as f:
		filedata = f.read()
	filedata = filedata.replace(',', ' ')
	with open(filename, "w") as f:
		f.write(filedata)	
if __name__ == '__main__':
	filename = "CORTEXM0DS_pads.VDD.pw_hl.pna"
	outputfile = "armcore.sp"
	raw_branch,all_layer, layer_info = read_line(filename)	
	branch_list = reform_branch(raw_branch, all_layer, layer_info)
	layer_name, branch_layer, layer_type = seperate_layer(branch_list)# 0 is layer 37, metal 9. 3 is layer 26, via 8, 4 is are wires
		
	tree_layer, tree_layername = get_tree(layer_name, branch_layer,layer_type) # output are trees on each metal layer. 
	nodes_position, node_current = collect_nodes(branch_list)
	write_spice(layer_name, branch_layer, tree_layer, layer_type, nodes_position, node_current,outputfile) 
	
