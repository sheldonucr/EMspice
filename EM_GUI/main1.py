#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 05:21:49 2019

@author: yibo
"""
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons
import numpy as np
import IPython
from IPython import get_ipython

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}

matplotlib.rc('font', **font)

import colorsys
color_cold=(0,0.0,0) # Select Coldest Color in RGB
color_hot=(0,0,0.0) # Select Hotest Color in HSV
(h1,s1,v1)=colorsys.rgb_to_hsv(color_cold[0],color_cold[1],color_cold[2])
(h2,s2,v2)=colorsys.rgb_to_hsv(color_hot[0],color_hot[1],color_hot[2])
def define_color(stress): #display the stress to color
    #apply HSV instead of RGB here, only HSV can change the color gradually. 
    #priting color of line segment will vary between the selected Cold and Hot color. 
    if (stress>500):
        stress=500
    if (stress<-500):
        stress=-500
    (r,g,b)=colorsys.hsv_to_rgb(h1+(h2-h1)*(stress+500)/1000,s1+(s2-s1)*(stress+500)/1000,v1+(v2-v1)*(stress+500)/1000)
    return((r,g,b))
def define_color2(current): #display the current density to color
    if (current>10000):
        current=10000
    if(current<-10000):
        current=-10000
    (r,g,b)=colorsys.hsv_to_rgb(h1+(h2-h1)*(current+10000)/20000,s1+(s2-s1)*(current+10000)/20000,v1+(v2-v1)*(current+10000)/20000)
    return((r,g,b))



def read_line(line_file):
    a = np.loadtxt(line_file).astype(np.float)
    return(a);
    
    
b=read_line("temp2/output.txt")
b2=read_line("temp2/output2.txt")
b3=read_line("temp2/stress.txt")
b4=read_line("temp2/void.txt")
b5=read_line("temp2/curden.txt")
m=b
m2=b2
m3=b3
m4=b4
m5=b5
line_number=np.size(m,0)
branch_number=np.size(m2,0)
stress_number=np.size(m3,0)
void_number=np.size(m4,0)
curden_number=np.size(m5,0)
match=0
for i in range(line_number):
    for j in range(branch_number): 
        if((m[i,4]==m2[j,3])and(m[i,5]==m2[j,4])and(m[i,6]==m2[j,5])and(m[i,7]==m2[j,6])):
            m[i,1]=m2[j,1]
            m[i,2]=m2[j,2]
for i in range(line_number):
    for j in range(stress_number):
        if((m[i,0]==m3[j,0])and(m[i,1]==m3[j,1])and(m[i,2]==m3[j,2])):
            m[i,11]=1
for i in range(line_number):
    for j in range(void_number):
        if((m[i,0]==m4[j,0])and(m[i,1]==m4[j,1])and(m[i,2]==m4[j,2])):
            m[i,13]=m4[j,5]
            match=match+1
            if((m[i,4]==m4[j,3])and(m[i,5]==m4[j,4])):
                m[i,12]=0
            else:
                m[i,12]=1
for i in range(line_number):
    for j in range(curden_number):
        if((m[i,0]==m5[j,0])and(m[i,1]==m5[j,1])and(m[i,2]==m5[j,2])):
            m[i,9]=m5[j,3]
            m[i,10]=1



    


#image alignment.
ipython = get_ipython()
ipython.magic("%matplotlib qt5") # inline or qt5
plt.rcParams ['figure.dpi']=100 #This is the default dpi 100.
plt.rcParams['figure.figsize']=(10.0,10.0) #default size is (6.0, 4.0), (6.0,6.0)will keep the correct length/width ratio

# count how many different layers
final_void=0
layer_label=[]
for i in range(line_number):
    if(not(m[i,0] in layer_label)):
       layer_label.append(m[i,0]) 

#plot the figure
fig, ax = plt.subplots()
l=[]
for i in range(line_number):
    layer_id=m[i,0]#m[i,0:13] is: layerID-treeID-BranchID-lineID-x1-y1-x2-y2-width-current_density-current_bit-stress_bit-void_bit - void_width
    tree_id=m[i,1]
    branch_id=m[i,2]
    line_id=m[i,3]#not necessary
    x1=m[i,4]
    y1=m[i,5]
    x2=m[i,6]
    y2=m[i,7]
    w=m[i,8]
    cd=m[i,9]
    c_bit=m[i,10]
    s_bit=m[i,11]
    v_bit=m[i,12]
    v_width=m[i,13]
    #if((c_bit==1)or(s_bit==1)or(v_bit==1)):
    
    '''line_stress = np.loadtxt(('stress'+str(l_id)+'.txt')).astype(np.int)'''# stress-[id].txt: an arrary, contanis the stress of ith line. 
    '''line_stress=[100, 200, -300, 400, 250, -100, 200]''' #Fake stress
    line_stress=[cd] #Fake stress
    num_segment = len(line_stress) #number of line segments    
    if((v_bit==-1)or((v_bit!=-1)and(v_width==0))):#no void or void length==0, just plot the branch without enlarging.
        if(x1==x2):
            for j in range(num_segment):
                x,=ax.fill([x1-w/2,x1+w/2,x2+w/2,x2-w/2],[y1+(y2-y1)/num_segment*j,y1+(y2-y1)/num_segment*j,y1+(y2-y1)/num_segment*(j+1),y1+(y2-y1)/num_segment*(j+1)],color=(0.6,0.6,0.6),visible=True,label=str(int(layer_id)))
                l.append(x) #x is the pointer to a ploygen
        else:
            for j in range(num_segment):
                x,=ax.fill([x1+(x2-x1)/num_segment*j,x1+(x2-x1)/num_segment*j,x1+(x2-x1)/num_segment*(j+1),x1+(x2-x1)/num_segment*(j+1)],[y1+w/2,y1-w/2,y2-w/2,y2+w/2],color=(0.6,0.6,0.6),visible=True,label=str(int(layer_id)))
                l.append(x)
    else:# the v_bit=0 or v_bit=1 means there is a void on the branch
        final_void=final_void+1
        w=w*40
        v_width=v_width*15
        #First plot the branch
        if(x1==x2):
            for j in range(num_segment):
                x,=ax.fill([x1-w/2,x1+w/2,x2+w/2,x2-w/2],[y1+(y2-y1)/num_segment*j,y1+(y2-y1)/num_segment*j,y1+(y2-y1)/num_segment*(j+1),y1+(y2-y1)/num_segment*(j+1)],color=(0.5,1,0.5),visible=True,label=str(int(layer_id)))
                l.append(x) #x is the pointer to a ploygen
        else:
            for j in range(num_segment):
                x,=ax.fill([x1+(x2-x1)/num_segment*j,x1+(x2-x1)/num_segment*j,x1+(x2-x1)/num_segment*(j+1),x1+(x2-x1)/num_segment*(j+1)],[y1+w/2,y1-w/2,y2-w/2,y2+w/2],color=(0.5,1,0.5),visible=True,label=str(int(layer_id)))
                l.append(x)
        #Then plot the Void
        if(x1==x2):
            if(v_bit==0):
                y2=y1+v_width
            else:
                y1=y2-v_width
            x,=ax.fill([x1-w/2,x1+w/2,x2+w/2,x2-w/2],[y1,y1,y2,y2],color=(0,0,0),visible=True,label=str(int(layer_id)))
            l.append(x)
        else:
            if(v_bit==0):
                x2=x1+v_width
            else:
                x1=x2-v_width
            x,=ax.fill([x1,x1,x2,x2],[y1-w/2,y1+w/2,y2+w/2,y2-w/2],color=(0,0,0),visible=True,label=str(int(layer_id)))
            l.append(x)


                    
# Make checkbuttons with all plotted lines with correct visibility
lines=l
rax = plt.axes([0.9, 0.35, 0.1, 0.3]) #the location of button bar [rightleft_x,rightleft_y,width,height]
labels = [str(line.get_label()) for line in lines]#the label of every single polygen.
visibility = [line.get_visible() for line in lines]
#check = CheckButtons(rax, labels, visibility)
# Attention, the label in the checkbutton's label=['25','26,'27''] must be same with the label of polygen, which in the labels=[str(line.get_label()) for line in lines]
'''
check = CheckButtons(rax, ['41','25','26','27'], [True,True,True,True]) 
''' # dynamically deside how many check buttons
button_label_list=[]
button_default_value=[]
for i in range(len(layer_label)):
    button_label_list.append(str(int(layer_label[i])))
    button_default_value.append(True)
check=CheckButtons(rax,button_label_list,button_default_value)




def func(label): #此处输入的label，是button上的名字string,用它来检索每条线的label
    '''
    index = labels.index(label)
    lines[index].set_visible(not lines[index].get_visible())
    plt.draw()'''
    #now, change index to index[], which means a group of indexes
    index=[]
    line_num=len(lines)
    for j in range(line_num):
        if(labels[j]==label):
            index.append(j)
    #print(index)
    for k in range(len(index)):
        lines[index[k]].set_visible(not lines[index[k]].get_visible())
    plt.draw()

check.on_clicked(func)
im=plt.show()


