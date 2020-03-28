#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#path='/home/yibo/EM GUI tool/raw_output.txt'

file1=open('data/chiptop.sp',encoding='utf-8')
file2=open('temp2/output2.txt','w',encoding='utf-8')


count=0
total_branch=0
for line in file1:
    count=count+1
    if(line[0]=='R'):# 'line' here is a string = # 1 3 3.57e+01 ....., need to convert it to list, strip by ' '.
        total_branch=total_branch+1
        line_list=line.split( ) # 'line_list' is a list of small strings=['R41_1_2', 'n1_1620161_481040', n1_1620161_480880, 2.8e-05]
        branch=line_list[0].split('-') #branch is a list of string=['R41','1','2']
        branch0=branch[0].split('R')#branch is a list of string=['','41']
        branch[0]=branch0[1]#now branch is a list of string=['41','1','2'], which is [layer_id, tree_id, branch_id]
        xy1=line_list[1].split('_')#xy1 is a list of string=['n1','1620161','481040']
        xy1.pop(0)
        xy2=line_list[2].split('_')#xy1 is a list of string=['n1','1620161','481040']
        xy2.pop(0)
        
        #new_line is a list=[ layer_id, tree_id,branch_id,x1, y1, x2, y2, current_density], total [0..7]
        new_line=[int(branch[0]), int(branch[1]),int(branch[2]),int(xy1[0])/1000,int(xy1[1])/1000,int(xy2[0])/1000,int(xy2[1])/1000,float(line_list[3])]
       

        for i in range(len(new_line)):
           file2.write(str(new_line[i]))
           file2.write(' ')
        file2.write('\n')
       
file1.close()
file2.close()
