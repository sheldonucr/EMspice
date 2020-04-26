# EMspice -- The full-chip level EM-IR coupled analysis for EM sign-off analysis

EMspice performs the mult-physics electrical and stress analysis of multi-segment interconnect wires. The repot consists of phython version and matlab version impementation of EMspice. 

**ICC to spice code in "EMspice_python/icc2spice/"**

**PG_solver code in EMspice_python/pg_sim/"**

**Coupled simulation code and framework in "EMspice_python/emspice_python/"**

**GUI code in "EMspice_python/em_gui/"**

Here are details on how to run the code

**Step 1. Get the spice file from ICC result**

This step we convert the output from ICC to spice file which can be read by PG_solver

Find the {pnafilename}.pna file from icc simulation. Run the syn2spice.py to get the spice file. 

Input is {pnafilename}.pna file

Output is {spicefilename}.sp file

Command is python syn2spice.py

**Step 2. Coupled simulation using EM slover and PG solver**
 
This step we do the coupled simulation with PG_solver and EM_solver

Use the spice file ({spicefilename}.sp) from the step 1, and the constant for simulation in the file EM_spice_constant.txt 
as the input for this step. EM_spice.py  is the EM solver and em_cmd is the binary code of PG solver. 
Do the filter after we get the current density from PG solver and then do coupled simulation with EM solver 
and PG solver. The stress, void and resistance are updated after each simulation step. 

Input is Spice({spicefilename}.sp) file(from step 1) EM_spice_constant.txt (constant for simulation).

Output is Void.txt: position of void. Lvoid.txt: Length of void. Length of void at corresponding position. 
u_stress.txt: stress on each segment tree_node_voltage.txt nodal voltage on all points.

Command is python EM_spice.py {spicefilename}.sp

**Step 3 Display the information in GUI**

This step display the simulation result from the GUI. 

The GUI part take structure information current density stress and void location and size information 
and display that to the user graphally. 

Input is {pnafilename}.pna file and {spicefilename}.sp file for geomotery information, u_stress for stress information, u_curden.txt for current density information, Void and LVoid for void information and tree_node_voltage.txt for voltage information. 

Output is displayed to user with MatplotLib. 

More information please refer to detailed guidlines in "EMspice_python/docs/"
