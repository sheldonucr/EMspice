#!/usr/bin/env python3
# -*- coding: utf-8 -*-

file1=open('data/u_curden_20.txt',encoding='utf-8')
file2=open('temp2/curden.txt','w',encoding='utf-8')


count=0
for line in file1:
    count=count+1
    if(line[0]=='R'):# 'line' here is a string

        line_list=line.split( ) # 'line_list' is a list of small strings=['R41_1_2', 'n1_1620161_481040', n1_1620161_480880, 2.8e-05]
        branch=line_list[0].split('-') #branch is a list of string=['R41','1','2']
        branch0=branch[0].split('R')#branch0 is a list of string=['','41']
        branch[0]=branch0[1]#now branch is a list of string=['41','1','2'], which is [layer_id, tree_id, branch_id]
        for i in range(3):
           file2.write(str(branch[i]))
           file2.write(' ')
        file2.write(str(int(float(line_list[1])/1000000)))
        file2.write('\n')
      

       

       
file1.close()
file2.close()

