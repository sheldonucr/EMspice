file1=open('data/ChipTop.pna',encoding='utf-8')
layer_list=[]
layer_count={}
layer_tyep_list={}
for line in file1:
    if(line[0]=='#'):# 'line' here is a string = # 1 3 3.57e+01 ....., need to convert it to list, strip by ' '.
       line_list=line.split( ) # 'line_list' is a list of small strings=['#', '1', '3.57e+01', ...]
       #new_line is a list=[ layer_id, tree_id,branch_id,line_id, x1, y1, x2, y2, width,current_value,current_bit,stress_bit,void_bit]
       x1=line_list[10]
       x2=line_list[12]
       y1=line_list[11]
       y2=line_list[13]
       layer=line_list[8]
       if (not(layer in layer_list)):
           layer_list.append(layer)
           if (x1==x2):
               layer_count[str(layer)+'H']=0
               layer_count[str(layer)+'V']=1
           else:
               layer_count[str(layer)+'H']=1
               layer_count[str(layer)+'V']=0
       else:
           if (x1==x2):
               layer_count[str(layer)+'V']=layer_count[str(layer)+'V']+1
           else:
               layer_count[str(layer)+'H']=layer_count[str(layer)+'H']+1
file1.close()
        

    

file1=open('data/ChipTop.pna',encoding='utf-8')
file2=open('data/ChipTop2.pna','w',encoding='utf-8')
for line in file1:
    if(line[0]=='#'):
       line_list=line.split( )
       x1=line_list[10]
       x2=line_list[12]
       y1=line_list[11]
       y2=line_list[13]
       layer=line_list[8]
       if (layer_count[str(layer)+'H']>layer_count[str(layer)+'V']):
           if(y1==y2):
               file2.write(line)
       if (layer_count[str(layer)+'H']<layer_count[str(layer)+'V']):
           if(x1==x2):
               file2.write(line)
    else:
        file2.write(line)

file1.close()
file2.close()
