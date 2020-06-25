import random

## some parameters
layer1 = 27
layer2 = 25
layer3 = 41
# num of nodes
row = 34				#row numbers
column = 32				#column numbers
vsource = 1.8
vsnum = 2
#iprint = 7e-04			#current source value
iprint = random.uniform(1,10)*10e-05			#current source value

## generate netlist file
fout = open('b.sp', 'w')
#fout.write('* ' + str(row) + ' x ' + str(column) + ' power network\n')
# voltage sources
for i in range(vsnum):
	fout.write('v%d n1_8000_%d 0 %g \n' % (i, i*5*72000+174000, vsource))
#fout.write('*\n')
# resistors
for i in range(column):
#	fout.write('** column: ' + str(i+1) + '\n')
	r_value2 = random.uniform(0.15,0.68)
	for j in range(row-1):
		lnode = 'n1' + '_' + str(i*48000+8000) + '_' + str(j*72000+174000)
		rnode = 'n1' + '_' + str(i*48000+8000) + '_' + str(j*72000+246000)
		rname = str(layer2) + '-' + str(i+100) + '-' + str(j+100)
		fout.write('R' + rname + ' ' + lnode + ' ' + rnode + ' ' + str(r_value2) + '\n')
for i in range(row):
#	fout.write('** row: ' + str(i+1) + '\n')
	r_value1 = random.uniform(0.15,0.68)
	for j in range(column-1):
		lnode = 'n1' + '_' + str(j*48000+8000) + '_' + str(i*72000+174000)
		rnode = 'n1' + '_' + str(j*48000+56000) + '_' + str(i*72000+174000)
		rname = str(layer1) + '-' + str(i+100) + '-' + str(j+100)
		if(i == 0):
			fout.write('R' + rname + ' ' + lnode + ' ' + rnode + ' ' + str(r_value1) + '\n')
		else:
			fout.write('R' + rname + ' ' + lnode + ' ' + rnode + ' ' + str(r_value1) + '\n')

# current sources
#fout.write('*\n')
inum = 0
for i in range(row):
	for j in range(column):
		fout.write('iB00_%d_v n1_%d_%d 0 %g \n' % (inum, j*48000+8000,i*72000+174000, iprint))
		inum += 1

# netlist file ends
fout.write('*\n.op\n')
fout.write('.print tran ')
for i in range(row):
	for j in range(column):
		fout.write('v(n%d_%d_%d) ' % (1, j*48000+8000, i*72000+174000))
		if(layer1 != layer2):
			fout.write('v(n%d_%d_%d) ' % (1, j*48000+8000, i*72000+174000))
#	if(i < vsnum):
#		fout.write('v(_X_n%d_8_%d) ' % (1, i*72+174))
fout.write('\n.end\n')
fout.close()

