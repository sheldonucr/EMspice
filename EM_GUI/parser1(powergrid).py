#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#path need to be added later!
#path='/home/yibo/EM GUI tool/raw_output.txt'

file1=open('data/ChipTop.pna',encoding='utf-8')
file2=open('temp2/output.txt','w',encoding='utf-8')

#line=file1.readline()
count=0
for line in file1:
    count=count+1
    if(line[0]=='#'):# 'line' here is a string = # 1 3 3.57e+01 ....., need to convert it to list, strip by ' '.
       line_list=line.split( ) # 'line_list' is a list of small strings=['#', '1', '3.57e+01', ...]
       #new_line is a list=[ layer_id, tree_id,branch_id,line_id, x1, y1, x2, y2, width,current_value,current_bit,stress_bit,void_bit]
       new_line=[int(line_list[8]), int(-1),int(-1),int(line_list[1]),int(line_list[10])/1000,int(line_list[11])/1000,int(line_list[12])/1000,int(line_list[13])/1000,float(line_list[9])/1000,float(line_list[5]),int(-1),int(-1),int(-1),int(-1)]
       

       for i in range(len(new_line)):
           file2.write(str(new_line[i]))
           file2.write(' ')
       file2.write('\n')
       
file1.close()
file2.close()
