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
from matplotlib.patches import Ellipse, Polygon
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 12}

matplotlib.rc('font', **font)



    
#Add a zoom-in/zoom-out function for mouse scroll 
def zoom_factory(ax,base_scale = 2.):
    def zoom_fun(event):
   
        cur_xlim = ax.get_xlim()# get the current x and y limits
        cur_ylim = ax.get_ylim()
        cur_xrange = (cur_xlim[1] - cur_xlim[0])*.5
        cur_yrange = (cur_ylim[1] - cur_ylim[0])*.5
        
        xdata = event.xdata # get event location_x and location_y
        ydata = event.ydata
        if event.button == 'up': # deal with zoom in
            scale_factor = 1/base_scale
        elif event.button == 'down': # deal with zoom out
            scale_factor = base_scale
        else:# deal with something that should never happen
            scale_factor = 1
            print(event.button)
        
        ax.set_xlim([xdata - cur_xrange*scale_factor, xdata + cur_xrange*scale_factor])# set new limits
        ax.set_ylim([ydata - cur_yrange*scale_factor, ydata + cur_yrange*scale_factor])
        plt.draw() # force re-draw

    fig = ax.get_figure() # get the figure of interest
    fig.canvas.mpl_connect('scroll_event',zoom_fun)   # attach the call back
    return zoom_fun#return the function


#Add a hot-key function for keyboard.

def key_factory(ax):
    def key_fun(event):
        print('keyboard function inside')
        print(event.key)
        cur_xlim = ax.get_xlim()# get the current x and y limits
        cur_ylim = ax.get_ylim()
        cur_xrange=(cur_xlim[1]-cur_xlim[0])
        cur_yrange=(cur_ylim[1]-cur_ylim[0])
        if event.key == "ctrl+c": # deal with zoom in
            print(' back to center')
            ax.set_xlim([xlim[0],xlim[1]])# set new limits
            ax.set_ylim([ylim[0], ylim[1]])
            plt.draw() # force re-draw
        if event.key=='up':
            ax.set_ylim([cur_ylim[0]+cur_yrange*0.1, cur_ylim[1]+cur_yrange*0.1]) 
            plt.draw() # force re-draw
        if event.key=='down':
            ax.set_ylim([cur_ylim[0]-cur_yrange*0.1, cur_ylim[1]-cur_yrange*0.1])
            plt.draw() # force re-draw
        if event.key=='left':
            ax.set_xlim([cur_xlim[0]-cur_xrange*0.1,cur_xlim[1]-cur_xrange*0.1])# set new limits
            plt.draw() # force re-draw
        if event.key=='right':
            ax.set_xlim([cur_xlim[0]+cur_xrange*0.1,cur_xlim[1]+cur_xrange*0.1])# set new limits
            plt.draw() # force re-draw
    fig = ax.get_figure() # get the figure of interest
    fig.canvas.mpl_connect('key_press_event',key_fun)   # attach the call back
    return key_fun#return the function



import colorsys
color_cold=(0,0,1) # Select Coldest Color in RGB
color_hot=(1,0,0) # Select Hotest Color in HSV
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
    (r,g,b)=colorsys.hsv_to_rgb(h1+(h2-h1)*(current)/25,s1+(s2-s1)*(current)/25,v1+(v2-v1)*(current)/25)
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
layer_label=[]
for i in range(line_number):
    if(not(m[i,0] in layer_label)):
       layer_label.append(m[i,0]) 
layer_label.sort()
total_layer=len(layer_label)

'''
#Add hatch pattern for layers.
pattern_lib=['/  ', '- ' , '+  ' , '  ' , '/-/-  ' ,'/\  ']
pattern_num=len(pattern_lib)
layer_pattern_dic={}
for i in range(total_layer):
    layer_pattern_dic.setdefault(layer_label[i],pattern_lib[i])
'''

#plot the figure
fig, ax = plt.subplots()
l=[]



final_match=0
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
    
    
    '''line_stress = np.loadtxt('temp/'+'stress'+str(layer_id)+'_'+str(tree_id)+'_'+str(branch_id)+'.txt').astype(np.int)'''
    #line_stress = np.loadtxt('stress'+str(layer_id)+'_'+str(tree_id)+'_'+str(branch_id)+'.txt').astype(np.int)
    '''line_stress=[100, 200, -300, 400, 250, -100, 200]''' #Fake stress

    if(not(s_bit==1)):#no stress warning, just plot the branch.
        line_stress=[0] #a fake stress value
        num_segment = len(line_stress) #number of line segments divided by stress
        if(x1==x2):
            for j in range(num_segment):
                #grid base color or stress color
                branch_color=(0.8,0.8,0.8)
                x,=ax.fill([x1-w/2,x1+w/2,x2+w/2,x2-w/2],[y1+(y2-y1)/num_segment*j,y1+(y2-y1)/num_segment*j,y1+(y2-y1)/num_segment*(j+1),y1+(y2-y1)/num_segment*(j+1)],color=branch_color,visible=True,label=str(int(layer_id)))
                l.append(x) #x is the pointer to a ploygen
                '''x=ax.add_patch(Polygon([[x1-w/2, y1], [x1+w/2, y1], [x2+w/2, y2], [x2-w/2, y2]], closed=True,fill=False, hatch=layer_pattern_dic[layer_id],label=str(int(layer_id))))
                l.append(x)'''
        else:
            for j in range(num_segment):
                branch_color=(0.8,0.8,0.8)
                x,=ax.fill([x1+(x2-x1)/num_segment*j,x1+(x2-x1)/num_segment*j,x1+(x2-x1)/num_segment*(j+1),x1+(x2-x1)/num_segment*(j+1)],[y1+w/2,y1-w/2,y2-w/2,y2+w/2],color=branch_color,visible=True,label=str(int(layer_id)))
                l.append(x)
                '''x=ax.add_patch(Polygon([[x1, y1+w/2], [x1, y1-w/2], [x2, y2-w/2], [x2, y2+w/2]], closed=True,fill=False, hatch=layer_pattern_dic[layer_id],label=str(int(layer_id))))
                l.append(x)'''   

    else:#Stress warning, load stress information
        pass;
special=[2819,
2821,
2835,
2836,
2855,
2857,
2871,
2872,
2891,
2893,
2907,
2908,
2932,
2933,
2952,
2954]
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
    
    
    '''line_stress = np.loadtxt('temp/'+'stress'+str(layer_id)+'_'+str(tree_id)+'_'+str(branch_id)+'.txt').astype(np.int)'''
    #line_stress = np.loadtxt('stress'+str(layer_id)+'_'+str(tree_id)+'_'+str(branch_id)+'.txt').astype(np.int)
    '''line_stress=[100, 200, -300, 400, 250, -100, 200]''' #Fake stress

    if(not(s_bit==1)):#no stress warning, just plot the branch.
        pass;
    else:#Stress warning, load stress information
        '''w=w*50'''
        line_stress = np.loadtxt('temp2/'+'stress'+str(int(layer_id))+'_'+str(int(tree_id))+'_'+str(int(branch_id))+'.txt').astype(np.int)
        num_segment = len(line_stress) #number of line segments divided by stress
        if(x1==x2):
            for j in range(num_segment):
                #grid base color or stress color
                branch_color=define_color(line_stress[j])
                x,=ax.fill([x1-w/2,x1+w/2,x2+w/2,x2-w/2],[y1+(y2-y1)/num_segment*j,y1+(y2-y1)/num_segment*j,y1+(y2-y1)/num_segment*(j+1),y1+(y2-y1)/num_segment*(j+1)],color=branch_color,visible=True,label=str(int(layer_id)))
                l.append(x) #x is the pointer to a ploygen
                '''x=ax.add_patch(Polygon([[x1-w/2, y1], [x1+w/2, y1], [x2+w/2, y2], [x2-w/2, y2]], closed=True,fill=False, hatch=layer_pattern_dic[layer_id],label=str(int(layer_id))))
                l.append(x)'''
        else:
            for j in range(num_segment):
                branch_color=define_color(line_stress[j])
                x,=ax.fill([x1+(x2-x1)/num_segment*j,x1+(x2-x1)/num_segment*j,x1+(x2-x1)/num_segment*(j+1),x1+(x2-x1)/num_segment*(j+1)],[y1+w/2,y1-w/2,y2-w/2,y2+w/2],color=branch_color,visible=True,label=str(int(layer_id)))
                l.append(x)
                '''x=ax.add_patch(Polygon([[x1, y1+w/2], [x1, y1-w/2], [x2, y2-w/2], [x2, y2+w/2]], closed=True,fill=False, hatch=layer_pattern_dic[layer_id],label=str(int(layer_id))))
                l.append(x)'''   
    '''if((v_bit!=-1)and(v_width!=0)):'''
    if(i in special):
        '''v_width=v_width*10'''
        v_width=1
        w=w*3
        final_match=final_match+1
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

#ax range data for the key_factory
xlim=ax.get_xlim()
ylim=ax.get_ylim()


#Adjust the following font(here means adjust the font of the button)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 10}

matplotlib.rc('font', **font)
# Make checkbuttons with all plotted lines with correct visibility
lines=l
rax = plt.axes([0.9, 0.35, 0.1, 0.3]) #the location of button bar [right_down_x,right_down_y,width,height]
labels = [str(line.get_label()) for line in lines]#the label of every single polygen.
visibility = [line.get_visible() for line in lines]

# Attention, the label in the checkbutton's label=['25','26,'27''] must be same with the label of polygen, which in the labels=[str(line.get_label()) for line in lines]
'''
check = CheckButtons(rax, ['41','25','26','27'], [True,True,True,True])#check = CheckButtons(rax, labels, visibility)
''' 
# dynamically deside how many check buttons
button_label_list=[]
button_default_value=[]
for i in range(len(layer_label)):
    #button_label_list.append(str(int(layer_label[i]))+' '+layer_pattern_dic[layer_label[i]])
    button_label_list.append('L'+str(int(layer_label[i])))
    button_default_value.append(True)
check=CheckButtons(rax,button_label_list,button_default_value)#check = CheckButtons(rax, labels, visibility)


'''
#Add corresponding graphic hatch example for checkbuttons
for i in range(total_layer):
    rax2 = plt.axes([0.952, 0.35+0.3/(total_layer+1)*(i+0.75), 0.047,0.3/(total_layer+1)*0.5 ]) #the location of button bar [right_down_x,right_down_y,width,height]
    rax2.set_axis_off() #Turn the x- and y-axis off.
    #rax2.fill([0,2,2,0],[0,0,1,1],hatch=pattern_lib[i])
    rax2.add_patch(Polygon([[0,0],[1,0],[1,1],[0,1]],closed=True,fill=False, hatch=layer_pattern_dic[layer_label[total_layer-1-i]]))
'''



def func(label): #此处输入的label，是button上的名字string,用它来检索每条线的label
    label1=label.split('L')
    label=label1[1]
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

#add colorbar
#Attention:  rax=plt.axes([bottomleft_x,bottomleft_y,width,height]) can not fully decide the colorbar size.
# The pyplot.colorbar() has a limited height range, which is 0.3 max.
# The pyplot.colorbar() has a fixed width, not adjustable.
rax = plt.axes([-0.182, 0.15, 0.2, 0.6])
sm = plt.cm.ScalarMappable(cmap="jet", norm=plt.Normalize(vmin=-500, vmax=500))
sm.set_array([])
plt.colorbar(sm)  

#Call the mouse zoom-in/out function.
scale=2
f=zoom_factory(ax,base_scale=scale)

#Call the keyboard function
f2=key_factory(ax)

#plot the image.
plt.show()

#
#input("Press <enter>")
