# EM-GUI
# 1. Files and Folders
## 1.1 Files
There are two types of files under the root folder:  six parser files and several function files.

![alt tag](https://github.com/zergzerg/EMspice_workrepo/blob/master/EM_GUI/EM_GUI%20matplotlib%20version/readme2.jpg "View of the files")

The parser files scan the original output information(such as the 'armcore.sp' ) and re-organize the power grid data.

Six Parser files:

     parser1(powergrid).py: Parse the  powergrid information from 'CORTEXM0DS_pads.VDD.pw_hl.pna' file.
     
     parser2(armcore).py: Parse the geometric information from 'armcore.sp'.
     
     parser3(stress).py: Parse the stress from 'u_stress_xx.txt'.
     
     parser4(curden).py: Parse the current density from 'u_curden_xx.txt'.
     
     parser5(void): Parse the void from 'u_void_xx.txt'.
     
     parser6(voltage).py: Parse the voltage from 'tree_node_voltage_xx.txt'
     

Each main.py file has its corresponding functions:

     main1.py:  Plot void.
     
     main2.py:  Plot stress.
     
     main3.py:  Plot current density.
     
     main4.py:  Plot stress+void.
     
     main5.py:  Plot current_denstiy+void.
     
     main6.py:  Plot stress+colorbar.
     
     main7.py:  Plot current_density+colorbar.
     
     main8.py:  Plot stress+void+colorbar.
    
     main9.py: Plot current_density+void+colorbar.
     
     main10.py: Plot void (const void_width).
     
     main11.py: Plot voltage.

## 1.2 Folders
 There are two folders: '/data'  and  '/temp2'.
 
 '/data' folder contains the original power grid data comes from the system.
 
 '/temp2' folder contains the temporary data of EM-GUI
 
 
 
 # 2. How To Run the EM-GUI
 
 ## Step1: Place power grid data in to  folder '/data'
Concretely, put 'CORTEXM0DS_pads.VDD.pw_hl.pna', 'armcore.sp', 'u_Lvoid_xx.txt', 'u_stress_xx.txt', 'u_curden_xx.txt', 'tree_node_voltage_xx' into the folder '/data'.

![alt tag](https://github.com/zergzerg/EMspice_workrepo/blob/master/EM_GUI/EM_GUI%20matplotlib%20version/readme3.jpg "no comment")
     
## (optional) Run 'H-V comperator'
This file eliminate some wrong branches in chiptop dataset. if you don't use chiptop, no need to do this.

## Step2: Run all parser files.
Run parser files one by one. Before executing each parser file, open the parser to make sure the input file name match the data file name which located in the '/data'.
Only run one parser between parser3 and parser6. They use same storage.  

![alt tag](https://github.com/zergzerg/EMspice_workrepo/blob/master/EM_GUI/EM_GUI%20matplotlib%20version/readme1.jpg "Adjust the path of input data")

## Step3: Use Ipython to run the main.py file.
  Select one of the main.py file, execute it with Ipython.

You only need to run parser files at the beginning, then you can execute different main.py files without running parser.py file again.

 
#For further help, contact the author: yliu401@ucr.edu / liuyibo1994@hotmail.com
 
 
 
 
