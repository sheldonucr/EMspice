# This is the EM solver coupled with pg solver from pg_sim

**Step 2. Coupled simulation using EM slover and PG solver**
 
This step we do the coupled simulation with PG_solver and EM_solver

Use the spice file ({spicefilename}.sp)  and the constant for simulation in the file EM_spice_constant.txt 
as the input for this step. EM_spice.py  is the EM solver and em_cmd is the binary code of PG solver. 
Do the filter after we get the current density from PG solver and then do coupled simulation with EM solver 
and PG solver. The stress, void and resistance are updated after each simulation step. 

Input is Spice({spicefilename}.sp) file(from step 1) EM_spice_constant.txt (constant for simulation).

Output is Void.txt: position of void. Lvoid.txt: Length of void. Length of void at corresponding position. 
u_stress.txt: stress on each segment tree_node_voltage.txt nodal voltage on all points.

Command is python EM_spice.py {spicefilename}.sp

